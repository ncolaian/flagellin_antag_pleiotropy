#! /usr/bin/env perl

# this script allows for the parrallelization of the MIT calculations

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
my $out = $ARGV[2];

#way to read all the contents of a file to a scalar
my $string = LoadFile("$out/sums.hash");
my %sums = %$string;
#^^^^^^^#### REALLY USEFUL ######^^^^^^^
open my $IN, "<", "$out/mit_mat_full.txt";

my $header = <$IN>;
chomp $header;
my @header_split = split /\t/, $header;

my %pos;
for ( my $i = 2; $i < scalar(@header_split); $i++ ) {
    $pos{$i} = split /-/, $header_split[$i];
}
open my $OUT, ">", "$out/calc_pt$start.txt";
my $line_num = 0;
while ( <$IN> ) {
    $line_num++;
    next if ( $line_num < $start || $line_num > $end );
    chomp $_;
    my @split_line = split /\t/, $_;
    #if the comparison is to itself skip MIT
    next if ( $split_line[0] == $split_line[1] );
    my $pos1 = shift @split_line;
    my $pos2 = shift @split_line;
    
    #get pos hrefs
    my $pos1_sum_href = $sums{$pos1};
    my $pos2_sum_href = $sums{$pos2};
    
    #now need to perform MIT
    my $total_comb = array_sum(\@split_line);
    
    my $total_mit = 0.000000;
    for ( my $i = 0; $i < scalar(@split_line); $i++ ) {
        next if ( $split_line[$i] == 0 );
        
        my @bases = $pos{$i+2};
        my $prob_1 = $pos1_sum_href->{$bases[0]}/$pos1_sum_href->{'total'};
        my $prob_2 = $pos2_sum_href->{$bases[1]}/$pos2_sum_href->{'total'};
        
        my $comb_prob = $split_line[$i]/$total_comb;
        
        $total_mit += $comb_prob * log($comb_prob/($prob_1*$prob_2));
    }
    print $OUT "$pos1\t$pos2\t$total_mit\n";
}
close($IN);
close($OUT);

#SUBROUTINE
sub array_sum {
    my ($aref) = @_;
    my $sum = 0;
    foreach my $part ( @$aref) {
        $sum += $part;
    }
    return $sum;
}