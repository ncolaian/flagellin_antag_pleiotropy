#! /usr/bin/env perl

# this script allows for the parrallelization of the matrix building

use strict;
use warnings;
use Readonly;
use Data::Dumper;
use YAML qw( DumpFile LoadFile );

Readonly::Array my @AMINO_ACIDS => qw/
A R N D C Q E G H I L K M F P S T W Y V
/;

Readonly::Array my @N_ACIDS => qw/ A T G C /;

my $start = $ARGV[0];
my $end = $ARGV[1];
my $type = $ARGV[2];
my $out = $ARGV[3];

#way to read all the contents of a file to a scalar
my $string = LoadFile("$out/build.hash");
my %hash = %$string;
##### REALLY USEFUL ######

my %out_pos;
if ( $type eq "a" ) {
    my $posit = 0;
    foreach my $aa ( @AMINO_ACIDS ) {
        foreach my $aa2 ( @AMINO_ACIDS ) {
            my $key = "$aa-$aa2";
            $out_pos{$key} = $posit;
            $posit++;
        }
    }
}

else {
    my $posit = 0;
    foreach my $na ( @N_ACIDS ) {
        foreach my $na2 ( @N_ACIDS ) {
            my $key = "$na-$na2";
            $out_pos{$key} = $posit;
            $posit++;
        }
    }
}

#Position by position analysis
open my $OUT, ">", "$out/mat_pt$start.txt";
for ( my $i = $start; $i <= $end; $i++ ) {
    for ( my $j = 1; $j <= keys %hash; $j++) {
        my @array;
        #I changed the 0 to 2 to handle low number of counts
        $array[(keys %out_pos)-1] = 0;
        @array = (0) x @array;
        my $pos = 0;
        foreach my $base ( @{$hash{$i}} ) {
            if ( $hash{$i}->[$pos] eq "-" || $hash{$j}->[$pos] eq "-" || $hash{$i}->[$pos] eq "X" || $hash{$j}->[$pos] eq "X" ) {
                $pos++;
                next;
            }
            #should put a try catch block here to handle errors *HANDLED ABOVE*
            my $tag = $hash{$i}->[$pos] . "-" . $hash{$j}->[$pos];
            $array[$out_pos{$tag}]++;
            $pos++;
        }
        print $OUT "$i\t$j\t";
        print $OUT join("\t", @array);
        print $OUT "\n";
    }
}

close($OUT);

#Small Subroutine
sub array_sum {
    my ($aref) = @_;
    my $sum = 0;
    foreach my $part ( @$aref) {
        $sum += $part;
    }
    return $sum;
}