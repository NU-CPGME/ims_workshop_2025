#!/usr/bin/perl

my $version = "0.2";

use strict;
use warnings;

my $usage = "
fastq_stats.pl  -  Output descriptive statistics on one or more fastq files.
		    Also works on gzipped files ending in .gz

usage: fastq_stats.pl [options] <reads_1.fastq(.gz)> <reads_2.fastq(.gz)> ...

options:
  -p	prefix. This option has no arguments. If this option is given, each
		argument, instead of pointing to an individual file, will instead
		represent a ([optional] path and) prefix of two paired read files and
		stats for both files combined will be output. The files must have
		the format:
			(/path/to/)<prefix>_1.fastq(.gz)
			(/path/to/)<prefix>_2.fastq(.gz)

";

use Getopt::Std;
our ($opt_p);
getopts('p');

die $usage unless (@ARGV);

my @results;
my @total;
for my $i (0 .. $#ARGV){
    my $skip;
	my @infiles;
	if ($opt_p){
		my $infile1 = "$ARGV[$i]\_1.fastq";
		my $infile2 = "$ARGV[$i]\_2.fastq";
		if (-e "$ARGV[$i]\_1.fastq.gz"){ #just presume that if one of the files is gzipped, the other will be too.
			$infile1 .= ".gz";
			$infile2 .= ".gz";
		}
		@infiles = ($infile1, $infile2);
	} else {
		@infiles = ($ARGV[$i]);
	}
	my $basename = $ARGV[$i];
	$basename =~ s/.*\///;
	$basename =~ s/_[12]\.fastq(.gz)*//;
	
	my @subresults;
	my @total;
	print STDERR "$basename\n";
	my $lineleng = 0;
    foreach my $infile (@infiles){
		my $in;
		if ($infile =~ m/\.gz$/){
			open ($in, "gunzip -cd $infile | ") or $skip = 1;
		} else {
			open ($in, "<$infile") or $skip = 1;
		}
		if ($skip){
			print STDERR "Could not open file $ARGV[$i]: $!\n";
			next;
		}
		my $readcount = 0;
		print STDERR "\e[K\r\tCounting read $readcount";
		my @lengths;
		my $first;
		while (<$in>){
			chomp(my $l1 = $_);
			if (!$first){
				next if ($l1 !~ m/^@/); #skip any header lines before the first read
				$first = 1;
			}
			chomp (my $l2 = <$in>);
			chomp (my $l3 = <$in>);
			chomp (my $l4 = <$in>);
			$readcount++;
			if ($readcount % 100000 == 0){
				print STDERR "\r\tCounting read $readcount";
			}
			my $leng = length($l2);
			push @lengths, $leng;
		push @total, $leng;
		}
		close ($in);
		my $line = "\r\tCounting read $readcount";
		$lineleng = length($line);
		print STDERR "$line";
		#my ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq) = stats(\@lengths);
		#push @results, ([$ARGV[$i], $sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq]);
	}
	my ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq) = stats(\@total);
	my @blanks = (" ") x $lineleng;
	print STDERR "\r", join("", @blanks);
	print STDERR "\r\t$num reads, $sum bases\n";
	push @results, ([$basename, $sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq]);
}
print STDERR "\n";
print "id\tnum_bases\tnum_reads\tmin_read_leng\tmax_read_leng\tmean_read_leng\tmedian_read_leng\tmode_read_leng\tmode_freq\n";
for my $i (0 .. $#results){
    my ($id, $sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq) = @{$results[$i]};
    print "$id\t$sum\t$num\t$min\t$maxi\t$rounded_mean\t$median\t$mode\t$mode_freq\n";
}

#--------------------
sub stats{
	my @lengths = @{$_[0]};
	my ($sum, $num, $min, $maxi, $mean, $median, $mode, $mode_freq);
	my %seen;
	my @sorted_leng = sort {$a <=> $b} @lengths;
	for my $i (0 .. $#sorted_leng){
		$sum += $sorted_leng[$i];
		$seen{$sorted_leng[$i]}++;
	}
	$num = $#sorted_leng + 1;
	$min = $sorted_leng[0];
	$maxi = $sorted_leng[$#sorted_leng];
	$mean = $sum/$num;
	my $rounded_mean = sprintf("%.2f", $mean);
	my @modes;
	foreach my $leng (sort {$seen{$b} <=> $seen{$a}} keys %seen) {
		push @modes, ([$leng, $seen{$leng}]);
	}
	$mode = $modes[0][0];
	$mode_freq = $modes[0][1];
	my $mid = int @sorted_leng/2;
	if (@sorted_leng % 2){
		$median = $sorted_leng[$mid];
	} else {
		$median = ($sorted_leng[$mid-1] + $sorted_leng[$mid])/2;
	}
	return ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq);
}

