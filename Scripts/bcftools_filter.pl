#!/usr/bin/perl

use strict;
use warnings;

my $usage = "

bcftools_filter.pl [options] <results.vcf or piped from bcftools>

Filters bcftools variant results produced from alignments against a reference.

Generate the input file by running:
samtools mpileup -E -M 0 -Q 25 -q 30 -m 2 -D -S -g -f <reference.fasta> <alignment.bam> | bcftools view -Icg - 
(for samtools and bcftools versions 0.1.19)

Output will be sequence file of reference strain where SNVs passing filters
are substituted into the sequence.  Missing data or filtered SNVs will be
replaced with gap characters '-'

Note: As a reminder, the default of samtools mpileup (without the '-B' flag) is
to perform BAQ base quality calculation. Though this can avoid calling false
SNPs around INDELs, it may result in some true bases matching the reference to
be filtered out of the output. Hence there may less false SNPs at the cost of
more false gaps. Caveat emptor.

Required:
  -f    fasta file for reference strain

Options:
  -q    minimum SNV quality score
        [default: 1]
  -c    minimum read consensus, in percent
        [default: 75]
  -d    minimum read depth
          * will be calculated from DP4 to only count reads above quality filter
        [default: 5]
  -D    maximum read depth, in fold of median for the isolate
        [default: 3]
  -r    minimum number of reads in each direction
        [default: 1]
  -h    DO NOT require homozygositiy
        [default: only keep SNVs homozygous under diploid model, i.e. GT = 1/1]
  -m    reference genome masking file, in interval format output by NCBI dustmaker,
        i.e.
        >chromosome_name
        209 - 215
        415 - 421
        487 - 494
        1104 - 1110
        etc.
  -o    Output sequence ID. If there is only one sequence in the reference, then
        its name will be replaced in the output. If there are two or more
        sequences in the reference, their IDs will be prefixed.
        (default: output will have reference sequence IDs)
        
";

use Getopt::Std;
use vars qw( $opt_f $opt_q $opt_c $opt_d $opt_D $opt_r $opt_h $opt_m $opt_o );
getopts('f:q:c:d:D:r:hm:o:');

die $usage unless @ARGV and $opt_f;

my $reffile     = $opt_f;
my $minqual     = $opt_q ? $opt_q : 200;
my $mincons     = $opt_c ? $opt_c : 75;
my $mindep      = $opt_d ? $opt_d : 5;
my $maxfold     = $opt_D ? $opt_D : 3;
my $mindir      = $opt_r ? $opt_r : 1;
my $maskfile    = $opt_m if $opt_m;
my $outid       = $opt_o if $opt_o;

my %mask;
if ($maskfile){
    my $totmask = 0;
    my $masklines = 0;
    open (my $in, "<$maskfile") or die "ERROR: Can't open $maskfile: $!\n";
    my $id;
    while (my $line = <$in>){
        chomp $line;
        if ($line =~ m/^>/){
            $id = substr($line, 1);
            $id =~ s/\s.*$//;
            next;
        }
        if ($line =~ m/^(\d+) - (\d+)/){
            $masklines++;
            my ($start, $stop) = ($1, $2);
            for my $i ($start .. $stop){
                $mask{$id}{$i} = 1;
                $totmask++;
            }
        }
    }
    close ($in);
    print STDERR "Total of $totmask masked positions in $masklines intervals.\n";
}

## Read in reference sequence
open (my $fin, "<$reffile") or die "ERROR: Can't open $reffile: $!\n";
my @seqs;
my %seqlengs;
my $totseqleng = 0;
my ($id, $seq);
while (my $line = <$fin>){
    chomp $line;
    next if $line =~ m/^\s*$/;
    if ($line =~ m/^>/){
        if ($id){
            push @seqs, ([$id, $seq]);
            my $leng = length($seq);
            $seqlengs{$id} = $leng;
            $totseqleng += $leng;
            print STDERR "$id: $leng bp\n";
            $seq = "";
        }
        $id = substr($line, 1);
        $id =~ s/\s.*$//;
        next;
    }
    $line =~ s/\s//g;
    $seq .= $line;
}
close ($fin);
if ($id){
    push @seqs, ([$id, $seq]);
    my $leng = length($seq);
    $seqlengs{$id} = $leng;
    $totseqleng += $leng;
    print STDERR "$id: $leng bp\n";
    $seq = "";
}
print STDERR "Total sequence length: $totseqleng\n";

#Read and filter the bcftools results
open (my $in, "<$ARGV[0]") or die "ERROR: Can't open $ARGV[0]: $!\n";
my %snvs; #contig, position, base, depth
my %missing_or_filtered;
#stats: below_min_qual(0), below_min_consensus(1), below_min_depth(2), above_max_depth(3), unidirectional(4), non-homozygous(5), masked(6), missing(7), filtered(8)
my @stats = (0) x 9;
my @depths;
my $last_id = " ";
my $last_pos = 0;
while (my $line = <$in>){
    chomp $line;
    next if $line =~ m/^\s*#/;
    my @tmp = split("\t", $line);
    my ($id, $pos, $alt, $qual, $stuff, $format, $format_val) = @tmp[0,1,4,5,7,8,9];
    if ($id ne $last_id){
        if ($last_pos != 0){
            my $last_leng = $seqlengs{$last_id};
            if ($last_pos < $last_leng){
                for my $i ($last_pos + 1 .. $last_leng){
                    $missing_or_filtered{$last_id}{$i} = 1;
                    $stats[7]++;
                    #push @depths, 0;
                }
            }
        }
        if ($pos > 1){
            for my $i (1 .. $pos - 1){
                $missing_or_filtered{$id}{$i} = 1;
                $stats[7]++;
                #push @depths, 0;
            }
        }
    }
    if ($last_pos > 0 and $last_pos + 1 != $pos){
        for my $i ($last_pos + 1 .. $pos - 1){
            $missing_or_filtered{$id}{$i} = 1;
            $stats[7]++;
            #push @depths, 0;
        }
    }
    my @formats = split(":", $format);
    my @format_vals = split(":", $format_val);
    my %fmt;
    for my $i (0 .. $#formats){
        $fmt{$formats[$i]} = $format_vals[$i];
    }
    $stuff =~ m/DP4=(\d+),(\d+),(\d+),(\d+)/;
    my ($rf, $rr, $af, $ar) = ($1, $2, $3, $4);
    my $total_depth;
    if ($rf or $rr or $af or $ar){
        $total_depth = $rf + $rr + $af + $ar;
        print STDERR "WARNING: DP4 does not match DP\n$line\n" if !$fmt{"DP"} or $fmt{"DP"} != $total_depth;
    } else {
        print STDERR "WARNING: No DP value in FORMAT section:\n$line\n" unless exists $fmt{"DP"};
        $total_depth = $fmt{"DP"};
    }
    push @depths, $total_depth unless $total_depth == 0;
    #push @depths, $total_depth;
    if ($alt ne "."){
        my $filt = 0;
        if ($qual < $minqual){
            $filt++;
            $stats[0]++;
        }
        my $cons = 100*(($af + $ar) / $total_depth);
        if ($cons < $mincons){
            $filt++;
            $stats[1]++;
        }
        if ($total_depth < $mindep){
            $filt++;
            $stats[2]++;
            $stats[7]++; #will also count these as missing / uncovered
        }
        if ($af < $mindir or $ar < $mindir){
            $filt++;
            $stats[4]++;
        }
        my $gtval;
        if ($fmt{"GT"}){
            $gtval = $fmt{"GT"};
        } else {
            print STDERR "WARNING: Variant position without GT value:\n$line\n";
        }
        if ($gtval ne "1/1" and $gtval ne "1" and !$opt_h){
            $filt++;
            $stats[5]++;
        }
        if ($mask{$id}{$pos}){
            $filt++;
            $stats[6]++;
        }
        if ($filt > 0){
            $missing_or_filtered{$id}{$pos} = 1;
            $stats[8]++;
        } elsif ($filt == 0) {
            @{$snvs{$id}{$pos}} = ($alt, $total_depth);
            if ($alt =~ m/,/){
                my $first = (split",", $alt)[0];
                @{$snvs{$id}{$pos}} = ($first, $total_depth);
            }
        }
    } else {
        #Eliminate non-SNP positions with high-qual coverage less than the minimum depth
        if ($total_depth < $mindep){
            $missing_or_filtered{$id}{$pos} = 1;
            $stats[7]++;
        }
    }
    $last_id = $id;
    $last_pos = $pos;
}
close ($in);
my $last_leng = $seqlengs{$last_id};
die "ERROR: No file given\n" unless $last_leng;
if ($last_pos < $last_leng){
    for my $i ($last_pos + 1 .. $last_leng){
        $missing_or_filtered{$last_id}{$i} = 1;
        $stats[7]++;
        #push @depths, 0;
    }
}

#calculate median depth
my $num_depths = scalar @depths;
print STDERR "num_depths = $num_depths\n";
@depths = sort{$a <=> $b} @depths;
my $mid = ($num_depths / 2) - 0.5;
my $median;
if ($mid == int($mid)){
    print STDERR "mid = $mid\n";
    $median = $depths[$mid];
} else {
    my $two = $num_depths / 2;
    my $one = $two - 1;
    $median = ($depths[$one] + $depths[$two]) / 2;
}
print STDERR "Median depth: $median\n";
print STDERR "Maximum depth: $depths[$#depths]\n";
my $maxdep = $median * $maxfold;
print STDERR "Maximum depth threshold: $maxdep\n";
#Filter SNVs above maximum depth
my $totalsnps = 0;
foreach my $id (keys %snvs){
    #print STDERR "$id\n";
    foreach my $pos (keys %{$snvs{$id}}){
        #print STDERR "^$pos\n";
        my ($alt, $depth) = @{$snvs{$id}{$pos}};
        #print STDERR "$id $pos $alt $depth\n";
        if ($depth > $maxdep){
            $stats[3]++;
            $stats[8]++;
            $missing_or_filtered{$id}{$pos} = 1;
            delete $snvs{$id}{$pos};
        } else {
            $totalsnps++;
        }
    }
}

print STDERR "Total missing or below minimum depth ($mindep): $stats[7]\n";
my $pct_covered = 100 * (($totseqleng - $stats[7]) / $totseqleng);
print STDERR "\tPercent aligned: $pct_covered\n";
print STDERR "Total filtered SNVs: $stats[8]\n";
#below_min_qual(0), below_min_consensus(1), below_min_depth(2), above_max_depth(3), unidirectional(4), non-homozygous(5), masked(6),
print STDERR "\tBelow minimum quality ($minqual): $stats[0]\n";
print STDERR "\tBelow minimum consensus ($mincons%): $stats[1]\n";
print STDERR "\tBelow minimum depth ($mindep): $stats[2]\n";
print STDERR "\tAbove maximum depth ($maxfold x $median): $stats[3]\n";
print STDERR "\tUnidirectional ($mindir): $stats[4]\n";
print STDERR "\tNon-homozygous: $stats[5]\n" if !$opt_h;
print STDERR "\tMasked: $stats[6]\n";
print STDERR "Total SNVs: $totalsnps\n";

print STDERR "\nOutputting replaced sequence\n";
foreach my $slice (@seqs){
    my ($id, $seq) = @{$slice};
    if ($outid){
        my $newid = "$outid\_$id";
        $newid = $outid if scalar @seqs == 1;
        print ">$newid\n";
    } else {
        print ">$id\n";
    }
    my @tmp = split(//, $seq);
    for my $i (0 .. $#tmp){
        my $base = $tmp[$i];
        if (my $snp = $snvs{$id}{$i+1}){
            $base = @{$snp}[0];
            #print STDERR "pos:". ($i + 1) . " SNV\n";
        }
        if ($missing_or_filtered{$id}{$i+1}){
            $base = "-";
            #print STDERR "pos:". ($i + 1) . " missing\n";
            #print STDERR "*** BOTH ***\n" if $snvs{$id}{$i+1};
        }
        print "$base";
    }
    print "\n";
}


