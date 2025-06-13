#!/usr/bin/perl

my $version = 0.2;

## Changes from v0.1
## Improve the efficiency of the contamination screen step

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use File::Spec::Functions qw ( catfile path );

#set defaults
my $cfile;
my $read1;
my $read2;
my $min_cov     = 5;
my $min_len     = 200;
my $bwa_loc;
my $blast_loc;
my $cont_file   = "phiX";
my $cont_id     = 98;
my $threads     = 8;
my $pref        = 'output';

my $usage = "
assembly_qc_filter.pl

1) Align reads back to assembly contigs to determine:
    a) average fold coverage per contig 
    b) minimum fold coverage per contig
    c) # of unaligned reads 
2) Remove contigs that align to phiX
    

Required:
  -c    assembly file of contigs in fasta format
  -1    forward read file
  -2    reverse read file
  
Optional:
  -m    Minimum average fold coverage of a contig to keep
        (default: $min_cov)
  -l    Minimum length of a contig to keep, in bases
        (default: $min_len)

  -a    path to bwa aligner executable
        (default: will search for bwa in PATH)
  -b    path to blastn executable
        (default: will search in PATH)
  -p    path to phiX sequence or other contaminant sequence file in fasta format
        (default: $cont_file)
  -x    minimum percent identity for phiX alingment
        (default: $cont_id)
  -t    threads
        (default: $threads)
  -o    output prefix
        (default: '$pref')

";

## Outputs:
## 1) fasta file of output sequences (minus contigs filtered by size or coverage)
## 2) statistics of input and output sequences (to STDOUT?)
##  number,total_size,average_size,median_size,minimum_size,maximum_size,n50,gc,average_coverage,median_coverage,maximum_coverage,minimum_coverage
## 3) per-contig statistics
##  id,size,gc,average_coverage,median_coverage,min_cov,max_cov

## command line processing.
use Getopt::Long;
Getopt::Long::Configure qw(gnu_getopt);

GetOptions(
    'c=s'           => \$cfile,
    '1=s'           => \$read1,
    '2=s'           => \$read2,
    'm=f'           => \$min_cov,
    'l=i'           => \$min_len,
    'a=s'           => \$bwa_loc,
    'b=s'           => \$blast_loc,
    'p=s'           => \$cont_file,
    'x=f'           => \$cont_id,
    't=i'           => \$threads,
    'o=s'           => \$pref
) or die "$usage";
die "$usage" unless $cfile and $read1 and $read2;

if ($bwa_loc){
    $bwa_loc = abs_path($bwa_loc);
    die "ERROR: bwa not found at $bwa_loc\n" unless (-x "$bwa_loc"); #check that it is executable
} else {
    #look in PATH
    my $path = is_path("bwa");
    if ($path){
        $bwa_loc = $path;
    } else {
        die "ERROR: bwa is required for optional read correction. No working version of bwa could be found in your PATH.\n";
    }
}

my %cont_key;
if ($cont_file =~ m/\S/){
    if ($blast_loc){
        $blast_loc = abs_path($blast_loc);
        die "ERROR: blastn not found at $blast_loc\n" unless (-x "$blast_loc"); #check that it is executable
    } else {
        #look in PATH
        my $path = is_path("blastn");
        if ($path){
            $blast_loc = $path;
        } else {
	    #### This is only for Quest
	    $blast_loc = "blastn";
            #die "ERROR: blastn is required for contamination filtering. No working version of blastn could be found in your PATH.\n";
        }
    }
    #this step is necessary for old blast (2.2.24)
    my $ccount = 0;
    if ($cont_file eq "phiX"){
        $cont_key{"Subject_0"} = "phiX";
    } else {
        open (my $cin, "<$cont_file") or die "ERROR: Can't open $cont_file: $!\n";
        while (my $line = <$cin>){
            chomp $line;
            if ($line =~ m/^>/){
                my $id = substr($line, 1);
                $id =~ s/\s.*$//;
                $ccount++;
                my $sid = "Subject_$ccount";
                $cont_key{$sid} = $id;
            }
        }
        close ($cin);
    }
}

## 1: count contigs in sequence file, get basic stats (num, total size, average size, median size, n50, gc, per-contig size)
#my %contigs;
print STDERR "Counting contigs...\n";
my @order;
my @contig_lengs;
my ($tot_gc, $tot_non_n) = (0) x 2;
my %starts;
my %stops;
open (my $cin, "<$cfile") or die "ERROR: Can't open $cfile: $!\n";
my $id;
my ($cleng, $cnon_n, $cgc) = (0) x 3;
my ($pre_num, $pre_size) = (0) x 2;
while (my $line = <$cin>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            @{$starts{$id}} = (0) x ($order[$#order][1] + 1);
            @{$stops{$id}} = (0) x ($order[$#order][1] + 1);
            push @contig_lengs, $order[$#order][1];
            $pre_size += $order[$#order][1];
            $tot_non_n += $order[$#order][2];
            $tot_gc += $order[$#order][3];
        }
        $id = substr($line, 1);
        $id =~ s/\s.*$//;
        push @order, ([$id, 0, 0, 0]);
        $pre_num++;
        #$contigs{$id}{'leng'} = 0;
        #$contigs{$id}{'non_n'} = 0;
        #$contigs{$id}{'gc'} = 0;
        next;
    }
    $line =~ s/\s//;
    $order[$#order][1] += length($line);
    $order[$#order][2] += ($line =~ tr/ACGTacgt/ACGTacgt/);
    $order[$#order][3] += ($line =~ tr/GCgc/GCgc/);
    $order[$#order][4] .= $line;
    
    #$contigs{$id}{'leng'} += length($line);
    #$contigs{$id}{'non_n'} += ($line =~ tr/ACGTacgt/ACGTacgt/);
    #$contigs{$id}{'gc'} += ($line =~ tr/GCgc/GCgc/);
}
if ($id){
    @{$starts{$id}} = (0) x ($order[$#order][1] + 1);
    @{$stops{$id}} = (0) x ($order[$#order][1] + 1);
    push @contig_lengs, $order[$#order][1];
    $pre_size += $order[$#order][1];
    $tot_non_n += $order[$#order][2];
    $tot_gc += $order[$#order][3];
} else {
    die "ERROR: contig file $cfile has no fasta records\n";
}
close ($cin);

@contig_lengs = sort{$a <=> $b}@contig_lengs;
my ($pre_min, $pre_max) = ($contig_lengs[0], $contig_lengs[$#contig_lengs]);
my ($pre_avg, $pre_stdev) = average(\@contig_lengs);
my $pre_med = median(\@contig_lengs);
my $pre_n50 = n50(\@contig_lengs);
my $pre_gcpct = sprintf("%.2f", 100 * ($tot_gc / $tot_non_n));

## 2: create bwa index from contigs
print STDERR "Indexing contigs...\n";
my $status = system("$bwa_loc index -p $pref\_tmp $cfile > /dev/null 2>&1");
die "ERROR: bwa index quit with status: $status\n" if $status;

## 3: align reads to contigs using bwa
## 4: count reads aligning per contig, store alignment starts and stops in separate hashes, then iterate over array to count starts and stop for coverage

open (my $bin, "$bwa_loc mem -t $threads $pref\_tmp $read1 $read2 2>/dev/null | ") or die "ERROR: Can't run bwa mem: $!\n";
my $unaligned = 0;
my $readcount = 0;
my $sec = 0;
my @tlens;
#open my $testout, ">test.sam";
print STDERR "\rAligning ($readcount)";
while (my $line = <$bin>){
    chomp $line;
    #print $testout "$line\n";
    next if $line =~ m/^@/;
    my @tmp = split("\t", $line);
    my ($flag, $id, $pos, $mapq, $cigar, $tlen) = @tmp[1..5,8];
    if ($flag & 256){
        $sec++;
        next;
        #skip secondary alignments
    }
    $readcount++;
    print STDERR "\rAligning ($readcount)" if $readcount % 1000 == 0;
    if ($flag & 4){
        $unaligned++;
        next;
    }
    if ($tlen > 0){
        #only use first template length
        push @tlens, $tlen;
    }
    ## filter for map quality?
    my @cigars = split(/(?<=\D)/, $cigar);
    my $end = $pos;
    foreach my $part (@cigars){
        if ($part =~ m/(\d+)([MD])/){
            my ($num, $op) = ($1, $2);
            $end += ($num - 1);
        }
    }
    next if $end == $pos;
    $starts{$id}[$pos]++;
    $stops{$id}[$end]++;
}
#close $testout;
close ($bin);
print STDERR "\rAligning ($readcount)\n";
my $pct_unaligned = sprintf("%.2f", 100 * ($unaligned / $readcount));
my ($avg_ins, $stdev_ins) = average(\@tlens);
print "#total_reads\tunaligned_reads\tpct_unaligned\tavg_ins_size\tstdev_ins_size\n";
print "$readcount\t$unaligned\t$pct_unaligned\t$avg_ins\t$stdev_ins\n\n";

#print STDERR "Total reads: $readcount, unaligned: $unaligned ($pct_unaligned%)\n";
#print STDERR "Average insert size: $avg_ins (stdev: $stdev_ins)\n\n";

unlink("$pref\_tmp.amb");
unlink("$pref\_tmp.ann");
unlink("$pref\_tmp.bwt");
unlink("$pref\_tmp.pac");
unlink("$pref\_tmp.sa");

my @results;
my @averages;

my @total_covs;
my $tout;
my %cont_blast_ids;
if ($cont_file =~ m/\S/){
    open ($tout, ">$pref.tmp.fasta");
}
print STDERR "Calculating coverages...\n";
foreach my $slice (@order){
    my ($id, $leng, $non_n, $gc, $seq) = @{$slice};
    my $gcpct = sprintf("%.2f", 100 * ($gc / $non_n));
    #print STDERR "$id (length $leng, gc $gcpct%)\n";
    my $cov = 0;
    my @covs;
    for my $i (1 .. $leng){
        $cov += $starts{$id}[$i];
        if ($cov < 0){
            #debug
            print STDERR "$id position:$i cov:$cov\n";
        }
        push @covs, $cov;
        push @total_covs, $cov;
        $cov -= $stops{$id}[$i];
    }
    @covs = sort{$a <=> $b}@covs;
    my ($min, $max) = ($covs[0], $covs[$#covs]);
    my ($avg, $stdev) = average(\@covs);
    my $med = median(\@covs);
    my @ads;
    foreach (@covs){
        push @ads, abs($_ - $med);
        #print STDERR "\t$_ ($ads[$#ads])\n" if $id eq "NODE_87_length_529_cov_26.0448";
    }
    my $mad = median(\@ads);
    push @results, ([$id, $leng, $gcpct, $avg, $stdev, $med, $mad, $min, $max, @covs]);
    push @averages, $avg;
    if ($leng >= $min_len and $avg >= $min_cov and $tout){
        my $outid = $id;
        $outid =~ s/\s.+//g;
        $cont_blast_ids{$outid} = $id;
        print $tout ">$outid\n$seq\n";
    }
    #print STDERR "\tavg:$avg median:$med ($mad) min:$min max:$max\n";
}
close $tout if $tout;
@total_covs = sort{$a <=> $b}@total_covs;
my ($pre_total_cov_min, $pre_total_cov_max) = ($total_covs[0], $total_covs[$#total_covs]);
my ($pre_total_cov_avg, $pre_total_cov_stdev) = average(\@total_covs);
my $pre_total_cov_med = median(\@total_covs);
#print STDERR "Average coverage: $pre_total_cov_avg (stdev:$pre_total_cov_stdev)\n";
#print STDERR "Median coverage: $pre_total_cov_med\n";

my %cont_filtered;
if ($cont_file =~ m/\S/ and -e "$pref.tmp.fasta"){
    print STDERR "Searching for contamination...\n";
    my $bin;
    if ($cont_file eq "phiX"){
        my $phixseq = phix();
        open ($bin, "echo \"$phixseq\" | $blast_loc -query $pref.tmp.fasta -subject - -outfmt 6 | ");
    } else {
        open ($bin, "$blast_loc -query $pref.tmp.fasta -subject $cont_file -outfmt 6 | ");
    }
    my $last_id = " ";
    my @over;
    while (my $line = <$bin>){
        chomp $line;
        my @tmp = split("\t", $line);
        my ($qid, $sid, $spct, $qstart, $qend, $sstart, $send) = @tmp[1,2,6,7,8,9];
        if ($qid ne $last_id){
            if (@over){
                my $full_id = $cont_blast_ids{$last_id};
                @over = sort{$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]}@over;
                my $outstring;
                foreach (@over){
                    my ($sid, $start, $stop, $pct) = @{$_};
                    if ($cont_key{$sid}){
                        $sid = $cont_key{$sid};
                    }
                    $outstring .= "$sid($start-$stop,$pct%) ";
                }
                chop $outstring;
                $cont_filtered{$full_id} = $outstring;
                @over = ();
            }
        }
        $last_id = $qid;
        unless ($spct < $cont_id) {
            push @over, ([$sid, $qstart, $qend, $spct]);
        }
    }
    close $bin;
    unlink "$pref.tmp.fasta";
    if (@over){
        my $full_id = $cont_blast_ids{$last_id};
        @over = sort{$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]}@over;
        my $outstring;
        foreach (@over){
            my ($sid, $start, $stop, $pct) = @{$_};
            if ($cont_key{$sid}){
                $sid = $cont_key{$sid};
            }
            $outstring .= "$sid($start-$stop,$pct%) ";
        }
        chop $outstring;
        $cont_filtered{$full_id} = $outstring;
        @over = ();
    }
}

#@averages = sort{$a <=> $b}@averages;
#my $med_avg = median(\@averages);
#my @ads;
#foreach (@averages){
#    push @ads, abs($_ - $med_avg);
#}
#my $mad_avg = median(\@ads);
#print STDERR "Median average: $med_avg (MAD:$mad_avg)\n";

## Output filtered sequences
print STDERR "Filtering and outputting ...\n";
open (my $out, ">$pref.filtered_sequences.fasta");
open (my $outc, ">$pref.contig_stats.txt");
print $outc "id\tlength\tgc_pct\tavg_cov\tstdev_cov\tmed_cov\tmin_cov\tmax_cov\tfiltered\n";
@total_covs = ();
@contig_lengs = ();
my @filtered_lengs;
my @filtered_covs;
my @contaminated_contigs;
my ($filt_size, $filt_num, $filt_gc, $filt_non_n) = (0) x 4;
($tot_gc, $tot_non_n) = (0) x 2;
my ($post_num, $post_size) = (0) x 2;
for my $i (0 .. $#results){
    my @tmp = @{$results[$i]};
    my ($id, $leng, $gcpct, $avg, $stdev, $med, $mad, $min, $max) = @tmp;
    print $outc "$id\t$leng\t$gcpct\t$avg\t$stdev\t$med\t$min\t$max\t";
    my @covs = @tmp[9..$#tmp];
    #print STDERR "$id ($leng, gc:$gcpct%) Average coverage: $avg ($min-$max)";
    my ($x1, $x2, $non_n, $gc, $seq) = @{$order[$i]};
    if ($leng >= $min_len and $avg >= $min_cov){
        ## check contig for contaminating sequence
        if ($cont_filtered{$id}){
            my $string = $cont_filtered{$id};
            print $outc "Y ($string)\n";
            push @contaminated_contigs, $id;
            next;
        }
        print $out ">$id\n$seq\n";
        $tot_gc += $gc;
        $tot_non_n += $non_n;
        push @contig_lengs, $leng;
        push @total_covs, @covs;
        $post_size += $leng;
        $post_num ++;
        print $outc "N";
    } else {
        $filt_gc += $gc;
        $filt_non_n += $non_n;
        push @filtered_lengs, $leng;
        push @filtered_covs, @covs;
        $filt_size += $leng;
        $filt_num ++;
        print $outc "Y";
    }
    print $outc "\n";
    #print STDERR "\n";
}
## print or output stats here
@contig_lengs = sort{$a <=> $b}@contig_lengs;
my ($post_min, $post_max) = ($contig_lengs[0], $contig_lengs[$#contig_lengs]);
my ($post_avg, $post_stdev) = average(\@contig_lengs);
my $post_med = median(\@contig_lengs);
my $post_n50 = n50(\@contig_lengs);
my $post_gcpct = sprintf("%.2f", 100 * ($tot_gc / $tot_non_n));
@total_covs = sort{$a <=> $b}@total_covs;
my ($post_total_cov_min, $post_total_cov_max) = ($total_covs[0], $total_covs[$#total_covs]);
my ($post_total_cov_avg, $post_total_cov_stdev) = average(\@total_covs);
my $post_total_cov_med = median(\@total_covs);

my ($filt_min, $filt_max, $filt_avg, $filt_stdev, $filt_med, $filt_n50, $filt_gcpct, $filt_total_cov_min, $filt_total_cov_max, $filt_total_cov_avg, $filt_total_cov_stdev, $filt_total_cov_med) = ("NA") x 12;
if (@filtered_lengs){
    @filtered_lengs = sort{$a <=> $b}@filtered_lengs;
    ($filt_min, $filt_max) = ($filtered_lengs[0], $filtered_lengs[$#filtered_lengs]);
    ($filt_avg, $filt_stdev) = average(\@filtered_lengs);
    $filt_med = median(\@filtered_lengs);
    $filt_n50 = n50(\@filtered_lengs);
    $filt_gcpct = sprintf("%.2f", 100 * ($filt_gc / $filt_non_n));
    @filtered_covs = sort{$a <=> $b}@filtered_covs;
    ($filt_total_cov_min, $filt_total_cov_max) = ($filtered_covs[0], $filtered_covs[$#filtered_covs]);
    ($filt_total_cov_avg, $filt_total_cov_stdev) = average(\@filtered_covs);
    $filt_total_cov_med = median(\@filtered_covs);
}

print STDERR "\n";

print "#pre/post\tnum\tsize\tavg\tstdev\tmedian\tmin\tmax\tn50\tgc\tavg_cov\tstdev_cov\tmedian_cov\tmin_cov\tmax_cov\n";
print "pre\t$pre_num\t$pre_size\t$pre_avg\t$pre_stdev\t$pre_med\t$pre_min\t$pre_max\t$pre_n50\t$pre_gcpct\t$pre_total_cov_avg\t$pre_total_cov_stdev\t$pre_total_cov_med\t$pre_total_cov_min\t$pre_total_cov_max\n";
print "post\t$post_num\t$post_size\t$post_avg\t$post_stdev\t$post_med\t$post_min\t$post_max\t$post_n50\t$post_gcpct\t$post_total_cov_avg\t$post_total_cov_stdev\t$post_total_cov_med\t$post_total_cov_min\t$post_total_cov_max\n";
print "filt\t$filt_num\t$filt_size\t$filt_avg\t$filt_stdev\t$filt_med\t$filt_min\t$filt_max\t$filt_n50\t$filt_gcpct\t$filt_total_cov_avg\t$filt_total_cov_stdev\t$filt_total_cov_med\t$filt_total_cov_min\t$filt_total_cov_max\n";
print "Contaminated contigs removed: ", scalar @contaminated_contigs, "\n";

#-------------------------------------------------------------------------------

sub is_path {
    ## Subroutine based on StackOverflow post by Sinan Unur (https://stackoverflow.com/a/8243770)
    my $exe = shift;
    my @path = path;
    my @pathext = ( q{} );
    if ($^O eq 'MSWin32'){
        push @pathext, map { lc } split /;/, $ENV{PATHEXT};
    }
    for my $dir (@path){
        for my $ext (@pathext){
            my $f = catfile $dir, "$exe$ext";
            return ($f) if -x $f;
        }
    }
    return();
}

sub average {
    my @array = @{$_[0]};
    my $num = scalar @array;
    return("NA","NA") unless $num > 0;
    my $sum = 0;
    foreach (@array){
        $sum += $_;
    }
    my $avg = $sum / $num;
    my $sos = 0;
    foreach (@array){
        $sos += ($_ - $avg)**2;
    }
    my $stdev = sqrt($sos / $num);
    return(sprintf("%.4f", $avg), sprintf("%.4f", $stdev));
    
}

sub median {
    my @array = @{$_[0]};
    my $leng = scalar @array;
    my $med;
    if ($leng % 2 == 0){
        my $mid2 = $leng / 2;
        my $mid1 = $mid2 - 1;
        $med = ($array[$mid1] + $array[$mid2]) / 2;
    } else {
        my $mid = ($leng / 2) - 0.5;
        $med = $array[$mid];
    }
    return($med);
}

sub n50 {
    my @array = @{$_[0]};
    my $size = 0;
    foreach (@array){
        $size += $_;
    }
    @array = sort {$b <=> $a} @array;
    my $count = 0;
    my $n50 = 0;
    foreach my $val (@array){
        $count += $val;
        if ($count >= $size/2){
            $n50 = $val;
            last;
        }
    }
    return($n50);
}

sub phix {
  my $phix_seq = ">NC_001422.1 Enterobacteria phage phiX174 sensu lato, complete genome
GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT
GATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAA
ATTCGACCTATCCTTGCGCAGCTCGAGAAGCTCTTACTTTGCGACCTTTCGCCATCAACTAACGATTCTG
TCAAAAACTGACGCGTTGGATGAGGAGAAGTGGCTTAATATGCTTGGCACGTTCGTCAAGGACTGGTTTA
GATATGAGTCACATTTTGTTCATGGTAGAGATTCTCTTGTTGACATTTTAAAAGAGCGTGGATTACTATC
TGAGTCCGATGCTGTTCAACCACTAATAGGTAAGAAATCATGAGTCAAGTTACTGAACAATCCGTACGTT
TCCAGACCGCTTTGGCCTCTATTAAGCTCATTCAGGCTTCTGCCGTTTTGGATTTAACCGAAGATGATTT
CGATTTTCTGACGAGTAACAAAGTTTGGATTGCTACTGACCGCTCTCGTGCTCGTCGCTGCGTTGAGGCT
TGCGTTTATGGTACGCTGGACTTTGTGGGATACCCTCGCTTTCCTGCTCCTGTTGAGTTTATTGCTGCCG
TCATTGCTTATTATGTTCATCCCGTCAACATTCAAACGGCCTGTCTCATCATGGAAGGCGCTGAATTTAC
GGAAAACATTATTAATGGCGTCGAGCGTCCGGTTAAAGCCGCTGAATTGTTCGCGTTTACCTTGCGTGTA
CGCGCAGGAAACACTGACGTTCTTACTGACGCAGAAGAAAACGTGCGTCAAAAATTACGTGCGGAAGGAG
TGATGTAATGTCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACT
AAAGGCAAGCGTAAAGGCGCTCGTCTTTGGTATGTAGGTGGTCAACAATTTTAATTGCAGGGGCTTCGGC
CCCTTACTTGAGGATAAATTATGTCTAATATTCAAACTGGCGCCGAGCGTATGCCGCATGACCTTTCCCA
TCTTGGCTTCCTTGCTGGTCAGATTGGTCGTCTTATTACCATTTCAACTACTCCGGTTATCGCTGGCGAC
TCCTTCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTA
CTGTAGACATTTTTACTTTTTATGTCCCTCATCGTCACGTTTATGGTGAACAGTGGATTAAGTTCATGAA
GGATGGTGTTAATGCCACTCCTCTCCCGACTGTTAACACTACTGGTTATATTGACCATGCCGCTTTTCTT
GGCACGATTAACCCTGATACCAATAAAATCCCTAAGCATTTGTTTCAGGGTTATTTGAATATCTATAACA
ACTATTTTAAAGCGCCGTGGATGCCTGACCGTACCGAGGCTAACCCTAATGAGCTTAATCAAGATGATGC
TCGTTATGGTTTCCGTTGCTGCCATCTCAAAAACATTTGGACTGCTCCGCTTCCTCCTGAGACTGAGCTT
TCTCGCCAAATGACGACTTCTACCACATCTATTGACATTATGGGTCTGCAAGCTGCTTATGCTAATTTGC
ATACTGACCAAGAACGTGATTACTTCATGCAGCGTTACCATGATGTTATTTCTTCATTTGGAGGTAAAAC
CTCTTATGACGCTGACAACCGTCCTTTACTTGTCATGCGCTCTAATCTCTGGGCATCTGGCTATGATGTT
GATGGAACTGACCAAACGTCGTTAGGCCAGTTTTCTGGTCGTGTTCAACAGACCTATAAACATTCTGTGC
CGCGTTTCTTTGTTCCTGAGCATGGCACTATGTTTACTCTTGCGCTTGTTCGTTTTCCGCCTACTGCGAC
TAAAGAGATTCAGTACCTTAACGCTAAAGGTGCTTTGACTTATACCGATATTGCTGGCGACCCTGTTTTG
TATGGCAACTTGCCGCCGCGTGAAATTTCTATGAAGGATGTTTTCCGTTCTGGTGATTCGTCTAAGAAGT
TTAAGATTGCTGAGGGTCAGTGGTATCGTTATGCGCCTTCGTATGTTTCTCCTGCTTATCACCTTCTTGA
AGGCTTCCCATTCATTCAGGAACCGCCTTCTGGTGATTTGCAAGAACGCGTACTTATTCGCCACCATGAT
TATGACCAGTGTTTCCAGTCCGTTCAGTTGTTGCAGTGGAATAGTCAGGTTAAATTTAATGTGACCGTTT
ATCGCAATCTGCCGACCACTCGCGATTCAATCATGACTTCGTGATAAAAGATTGAGTGTGAGGTTATAAC
GCCGAAGCGGTAAAAATTTTAATTTTTGCCGCTGAGGGGTTGACCAAGCGAAGCGCGGTAGGTTTTCTGC
TTAGGAGTTTAATCATGTTTCAGACTTTTATTTCTCGCCATAATTCAAACTTTTTTTCTGATAAGCTGGT
TCTCACTTCTGTTACTCCAGCTTCTTCGGCACCTGTTTTACAGACACCTAAAGCTACATCGTCAACGTTA
TATTTTGATAGTTTGACGGTTAATGCTGGTAATGGTGGTTTTCTTCATTGCATTCAGATGGATACATCTG
TCAACGCCGCTAATCAGGTTGTTTCTGTTGGTGCTGATATTGCTTTTGATGCCGACCCTAAATTTTTTGC
CTGTTTGGTTCGCTTTGAGTCTTCTTCGGTTCCGACTACCCTCCCGACTGCCTATGATGTTTATCCTTTG
AATGGTCGCCATGATGGTGGTTATTATACCGTCAAGGACTGTGTGACTATTGACGTCCTTCCCCGTACGC
CGGGCAATAACGTTTATGTTGGTTTCATGGTTTGGTCTAACTTTACCGCTACTAAATGCCGCGGATTGGT
TTCGCTGAATCAGGTTATTAAAGAGATTATTTGTCTCCAGCCACTTAAGTGAGGTGATTTATGTTTGGTG
CTATTGCTGGCGGTATTGCTTCTGCTCTTGCTGGTGGCGCCATGTCTAAATTGTTTGGAGGCGGTCAAAA
AGCCGCCTCCGGTGGCATTCAAGGTGATGTGCTTGCTACCGATAACAATACTGTAGGCATGGGTGATGCT
GGTATTAAATCTGCCATTCAAGGCTCTAATGTTCCTAACCCTGATGAGGCCGCCCCTAGTTTTGTTTCTG
GTGCTATGGCTAAAGCTGGTAAAGGACTTCTTGAAGGTACGTTGCAGGCTGGCACTTCTGCCGTTTCTGA
TAAGTTGCTTGATTTGGTTGGACTTGGTGGCAAGTCTGCCGCTGATAAAGGAAAGGATACTCGTGATTAT
CTTGCTGCTGCATTTCCTGAGCTTAATGCTTGGGAGCGTGCTGGTGCTGATGCTTCCTCTGCTGGTATGG
TTGACGCCGGATTTGAGAATCAAAAAGAGCTTACTAAAATGCAACTGGACAATCAGAAAGAGATTGCCGA
GATGCAAAATGAGACTCAAAAAGAGATTGCTGGCATTCAGTCGGCGACTTCACGCCAGAATACGAAAGAC
CAGGTATATGCACAAAATGAGATGCTTGCTTATCAACAGAAGGAGTCTACTGCTCGCGTTGCGTCTATTA
TGGAAAACACCAATCTTTCCAAGCAACAGCAGGTTTCCGAGATTATGCGCCAAATGCTTACTCAAGCTCA
AACGGCTGGTCAGTATTTTACCAATGACCAAATCAAAGAAATGACTCGCAAGGTTAGTGCTGAGGTTGAC
TTAGTTCATCAGCAAACGCAGAATCAGCGGTATGGCTCTTCTCATATTGGCGCTACTGCAAAGGATATTT
CTAATGTCGTCACTGATGCTGCTTCTGGTGTGGTTGATATTTTTCATGGTATTGATAAAGCTGTTGCCGA
TACTTGGAACAATTTCTGGAAAGACGGTAAAGCTGATGGTATTGGCTCTAATTTGTCTAGGAAATAACCG
TCAGGATTGACACCCTCCCAATTGTATGTTTTCATGCCTCCAAATCTTGGAGGCTTTTTTATGGTTCGTT
CTTATTACCCTTCTGAATGTCACGCTGATTATTTTGACTTTGAGCGTATCGAGGCTCTTAAACCTGCTAT
TGAGGCTTGTGGCATTTCTACTCTTTCTCAATCCCCAATGCTTGGCTTCCATAAGCAGATGGATAACCGC
ATCAAGCTCTTGGAAGAGATTCTGTCTTTTCGTATGCAGGGCGTTGAGTTCGATAATGGTGATATGTATG
TTGACGGCCATAAGGCTGCTTCTGACGTTCGTGATGAGTTTGTATCTGTTACTGAGAAGTTAATGGATGA
ATTGGCACAATGCTACAATGTGCTCCCCCAACTTGATATTAATAACACTATAGACCACCGCCCCGAAGGG
GACGAAAAATGGTTTTTAGAGAACGAGAAGACGGTTACGCAGTTTTGCCGCAAGCTGGCTGCTGAACGCC
CTCTTAAGGATATTCGCGATGAGTATAATTACCCCAAAAAGAAAGGTATTAAGGATGAGTGTTCAAGATT
GCTGGAGGCCTCCACTATGAAATCGCGTAGAGGCTTTGCTATTCAGCGTTTGATGAATGCAATGCGACAG
GCTCATGCTGATGGTTGGTTTATCGTTTTTGACACTCTCACGTTGGCTGACGACCGATTAGAGGCGTTTT
ATGATAATCCCAATGCTTTGCGTGACTATTTTCGTGATATTGGTCGTATGGTTCTTGCTGCCGAGGGTCG
CAAGGCTAATGATTCACACGCCGACTGCTATCAGTATTTTTGTGTGCCTGAGTATGGTACAGCTAATGGC
CGTCTTCATTTCCATGCGGTGCACTTTATGCGGACACTTCCTACAGGTAGCGTTGACCCTAATTTTGGTC
GTCGGGTACGCAATCGCCGCCAGTTAAATAGCTTGCAAAATACGTGGCCTTATGGTTACAGTATGCCCAT
CGCAGTTCGCTACACGCAGGACGCTTTTTCACGTTCTGGTTGGTTGTGGCCTGTTGATGCTAAAGGTGAG
CCGCTTAAAGCTACCAGTTATATGGCTGTTGGTTTCTATGTGGCTAAATACGTTAACAAAAAGTCAGATA
TGGACCTTGCTGCTAAAGGTCTAGGAGCTAAAGAATGGAACAACTCACTAAAAACCAAGCTGTCGCTACT
TCCCAAGAAGCTGTTCAGAATCAGAATGAGCCGCAACTTCGGGATGAAAATGCTCACAATGACAAATCTG
TCCACGGAGTGCTTAATCCAACTTACCAAGCTGGGTTACGACGCGACGCCGTTCAACCAGATATTGAAGC
AGAACGCAAAAAGAGAGATGAGATTGAGGCTGGGAAAAGTTACTGTAGCCGACGTTTTGGCGGCGCAACC
TGTGACGACAAATCTGCTCAAATTTATGCGCGCTTCGATAAAAATGATTGGCGTATCCAACCTGCA";
  return($phix_seq)
}