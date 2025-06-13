my $version = "0.1";

use strict;
use warnings;
use File::Basename;

$|++;

my $usage = "
fasta_list.pl  -   Lists sizes of all records in a fasta file

usage: fasta_list.pl [options] <fasta_file>

options:
  -s    Order output by size, largest to smallest
        [default: output order will be the same as the input]

";

use Getopt::Std;
our ($opt_s);
getopts('s');
die $usage unless (@ARGV);

open (my $in, "<$ARGV[0]") or die "ERROR: Can't open file $ARGV[0]: $!\n";
my @results;
my ($id, $seq);
while (my $line = <$in>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            my ($leng, $gc, $num_gap) = (0, "NA", "NA");
            if ($seq){
                $leng = length($seq);
                my $num_gc = $seq =~ tr/GgCc/GgCc/;
                my $num_at = $seq =~ tr/AaTt/AaTt/;
                $num_gap = $seq =~ tr/-Nn*/-Nn*/;
                my $base_tot = $num_gc + $num_at;
                $gc = sprintf("%.1f", 100 * ($num_gc/$base_tot));
            }
            push @results, ([$id, $leng, $gc, $num_gap]);
        }
        $id = substr($line, 1);
        $seq = "";
        next;
    }
    $line =~ s/\s//g;
    $seq .= $line;
}
close ($in);
if ($id){
    my ($leng, $gc, $num_gap) = (0, "NA", "NA");
    if ($seq){
        $leng = length($seq);
        my $num_gc = $seq =~ tr/GgCc/GgCc/;
        my $num_at = $seq =~ tr/AaTt/AaTt/;
        $num_gap = $seq =~ tr/-Nn*/-Nn*/;
        my $base_tot = $num_gc + $num_at;
        $gc = sprintf("%.1f", 100 * ($num_gc/$base_tot));
    }
    push @results, ([$id, $leng, $gc, $num_gap]);
}
die "ERROR: No records found in file\n" unless @results;
if ($opt_s){
    @results = sort{$b->[1] <=> $a->[1]} @results;
}
print "id\tlength\tpct_gc\tnum_gap\n";
foreach my $slice (@results){
    print join("\t", @{$slice}), "\n";
}
