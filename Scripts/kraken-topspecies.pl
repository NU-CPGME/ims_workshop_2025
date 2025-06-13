#!/usr/bin/perl

use strict;
use warnings;

my $usage = "
kraken-topspecies.pl <kraken-report results> <strain id>

Will output the species with the most read support

";

die $usage unless @ARGV;
my $strain = $ARGV[1] ? $ARGV[1] : "strain";

open (my $in, "<$ARGV[0]") or die "ERROR: Can't open $ARGV[0]: $!\n";
while (my $line = <$in>){
    chomp $line;
    if ($line =~ m/^\s+(\d+\.\d+)\s+\d+\s+\d+\s+(\S)\s+\d+\s+(.*)/){
        my ($pct, $type, $name) = ($1, $2, $3);
        if ($type eq "S"){
            print "$strain\t$name ($pct)\n";
            last;
        }
    }
}
close ($in);

