#!/opt/local/bin/perl

my $version = "0.2.2";

## Changes from v0.2.1
## Input genome size in bp rather than Mbp

## Changes from v 0.2
## Add option to skip gzipping output

## Changes from v 0.1:
##  - Uses getopts
##  - Option to enter expected genome size and desired average coverage. Script will calculate how many reads to output based on read lengths
##  - Added option to force re-output of all reads if less than the desired number are present in the input files.

use strict;
use warnings;

## Usage
my $usage = "
downsample_paired.pl  -  Randomly downsamples paired fastq files.  Outputs
                         gzipped fastq files.

usage:  downsample_paired.pl \\
        [options] \\
        <read file 1 [fastq or fastq.gz]> \\
        <read file 2 [fastq or fastq.gz]>

options:
  -p    Output prefix (default \"downsampled\");
  -g    Estimated genome size, in bp
  -f    Desired average coverage, in fold (i.e. 100 for 100-fold coverage)
  -l    Average length of reads in file1 and file2, separated by a comma (i.e. 101,99)
        If no values given, the average read length will be calculated (takes longer)
        If read lengths are given, be aware that the script does not check whether they are accurate.
  -r    Number of reads to output per file (i.e. 1000 will output 1000 R1 reads and 1000 R2 reads)
        If -r is given, -g, -f, and -l will be ignored.
        If -r is not given, values must be entered for -g and -f
  -o    Force output of all reads into new files if number of reads in original files is less than
        number of reads desired. (default is to not output any reads in this case)
  -z    Do not gzip output files (default is to gzip the output)

";
## command line processing.
use Getopt::Std;
use vars qw( $opt_p $opt_g $opt_f $opt_l $opt_r $opt_o $opt_z);
getopts('p:g:f:l:r:oz');
die "*****A value must be given either for -r or for both -g and -f\n\n$usage" if (!$opt_r and (!$opt_g or !$opt_f));
my $pref    = $opt_p ? $opt_p : "downsampled";
my $gensize = $opt_g if $opt_g;
my $fold    = $opt_f if $opt_f;
my $rleng   = $opt_l if $opt_l;
my $numseq  = $opt_r if $opt_r;

die $usage unless (@ARGV == 2);

my ($infile1, $infile2) = @ARGV;

die "****Two lengths separated by a comma must be given to -l\n\n$usage" if ($rleng and $rleng !~ m/\d+,\d+/);

#count lines
print STDERR "Counting reads\n";
my ($in, $line);
if ($infile1 =~ m/\.gz$/){
    if (!$opt_r and !$opt_l){
        open ($in, "gunzip -cd $infile1 | ");
    } else {
        chomp ($line = `gunzip -cd $infile1 | wc -l`);
    }
} else {
    if (!$opt_r and !$opt_l){
        open ($in, "<", $infile1);
    } else {
        chomp ($line = `wc -l $infile1`);
    }
}
my $count1;
if ($in){
    print STDERR "\tCalculating read lengths for file 1 ...";
    my @lengs;
    while (<$in>){
        my $id = $_;
        chomp(my $seq = <$in>);
        push @lengs, length($seq);
        my $id2= <$in>;
        my $qal = <$in>;
    }
    close ($in);
    ($count1, $rleng) = get_lengths(\@lengs);
    my $oleng = sprintf("%.2f", $rleng);
    print STDERR " avg read length: $oleng\n";
} else {
    $line =~ m/(\d+)/;
    $count1 = $1 / 4;
}
if ($infile2 =~ m/\.gz$/){
    if (!$opt_r and !$opt_l){
        open ($in, "gunzip -cd $infile2 | ");
    } else {
        chomp ($line = `gunzip -cd $infile2 | wc -l`);
    }
} else {
    if (!$opt_r and !$opt_l){
        open ($in, "<", $infile2);
    } else {
        chomp ($line = `wc -l $infile2`);
    }
}
my $count2;
if ($in){
    print STDERR "\tCalculating read lengths for file 2 ...";
    my @lengs;
    while (<$in>){
        my $id = $_;
        chomp(my $seq = <$in>);
        push @lengs, length($seq);
        my $id2= <$in>;
        my $qal = <$in>;
    }
    close ($in);
    my ($count, $avg) = get_lengths(\@lengs);
    $count2 = $count;
    $rleng .= ",$avg";
    my $oleng = sprintf("%.2f", $avg);
    print STDERR " avg read length: $oleng\n";
} else {
    $line =~ m/(\d+)/;
    $count2 = $1 / 4;
}
die "Read counts do not match!\n($infile1: $count1 reads, $infile2: $count2 reads)\n" if ($count1 != $count2);

#calculate $numseq if not given
if (!$numseq){
    my ($leng1, $leng2) = split(',', $rleng);
    $numseq = int(($gensize * $fold)/($leng1 + $leng2));
    print STDERR "\tNumber of reads per file needed to achieve ~".$fold."x coverage of $gensize bp genome: $numseq\n";
}

if ($count1 < $numseq){
    if ($opt_o){
        print STDERR "You only have $count1 reads in $infile1! All reads will be output\n";
        $numseq = $count1;
    } else {
        die "You only have $count1 reads in $infile1! Quitting.\n"
    }
}

#load array with read counts and randomize;
print STDERR "Randomizing\n";
my @read_nums;
for my $i (1 .. $count1){
    push @read_nums, $i;
}
shuffle_array(\@read_nums);
my @out_array;
for my $i (0 .. ($numseq - 1)){
    push @out_array, $read_nums[$i];
}
@out_array = sort {$a <=> $b} @out_array;

#output sequences
my ($in1, $in2);
if ($infile1 =~ m/\.gz$/){
    open ($in1, "gunzip -cd $infile1 | ") or die "Can't open input file 1: $!\n";
} else {
    open ($in1, "<", $infile1) or die "Can't open input file 1: $!\n";
}
if ($infile2 =~ m/\.gz$/){
    open ($in2, "gunzip -cd $infile2 | ") or die "Can't open input file 2: $!\n";
} else {
    open ($in2, "<", $infile2) or die "Can't open input file 2: $!\n";
}
my $index = 0;
my ($outstring1, $outstring2) = (" | gzip -c - > $pref\_1.fastq.gz", " | gzip -c - > $pref\_2.fastq.gz");
($outstring1, $outstring2) = ("> $pref\_1.fastq", "> $pref\_2.fastq") if $opt_z;
open (my $out1, $outstring1) or die "Can't open output file 1: $!\n";
open (my $out2, $outstring2) or die "Can't open output file 2: $!\n";
print STDERR "\rOutputting reads ";
for my $i (0 .. $#read_nums){
    my ($curRead1, $curRead2);
    for my $j (1 .. 4){
        $curRead1 .= <$in1>;
        $curRead2 .= <$in2>;
    }
    if ($i == $out_array[$index]){
        print $out1 "$curRead1";
        print $out2 "$curRead2";
        $index++;
    }
    print STDERR "\rOutputting reads $index" if $index % 10000 == 0;
    last if ($index == $numseq);
}
print STDERR "\rOutputting reads $index\n";
close ($in1);
close ($in2);
close ($out1);
close ($out2);

#-------------------------
sub get_lengths {
    my @array = @{$_[0]};
    my $count = @array;
    my $tot;
    for my $i (0 .. $#array){
        $tot += $array[$i];
    }
    my $avg = $tot / $count;
    return ($count, $avg);
}

sub shuffle_array {
    my $array = shift;
    return unless @$array;
    my $i = @$array;
    while (--$i){
        my $j = int rand ($i+1);
        @$array[$i,$j] = @$array[$j, $i];
    }
}
