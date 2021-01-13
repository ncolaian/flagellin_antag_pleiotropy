#! /usr/bin/env perl

# This script will get the location of each column in the full protein. I want to make sure that the positioning is comparable to all the files,
# so i am using flg22 as an anchor. All positioning is proportional to pos 13 in flg22 (the position that is D most often)

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Data::Dumper;

#My Variables
my $help = 0;
my $man = 0;
my $og_fasta;
my $trim_fasta;
my $out_file;

#Read in the variables from the command line
GetOptions( 'man|m'   =>  \$man,
            'help|h'  =>  \$help,
            'original_fasta|og=s'   => \$og_fasta,
            'trimal_fasta|tf=s'  => \$trim_fasta,
            'out_file|o=s' => \$out_file,
            ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manual and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

## Main ##

#In theory, I only need to perform the identification of the columns one at a time.
open my $IN, "<", $trim_fasta;
open my $OUT, ">", $out_file;

#get the first line without a -
my $id_line;
my $seq;
while ( my $line =  <$IN> ) {
    $seq = <$IN>;
if ( $seq !~ /-/ ) {
    $id_line = $line;
    last;
}
}
chomp $id_line;
$id_line = (split / /, $id_line)[0];
chomp $seq;

close $IN;

#print "$id_line\n$seq\n";
my $positioning_aref = get_positioning_of_cols($id_line, $seq);

print $OUT join("\n", @$positioning_aref);
close $OUT;

## Subroutines ##
sub get_positioning_of_cols {
   my ($id, $seq) = @_;
   
   my $orig = `grep -A 1 "$id" $og_fasta`;
   chomp $orig;
   print "$seq\n$orig\n";
   
   my $orig_seq = (split /\n/, $orig)[1];
   
   my @oseq_split = split //, $orig_seq;
   
   my @positions;
   my $spot = 0;
   foreach my $base ( split //, $seq ) {
      $spot = get_next_pos($spot,$base,@oseq_split);
      $spot = $spot + 1; #This gets next location in seq and actual position number
      push @positions, $spot;
   }
#   print Dumper(@positions);
   return(\@positions);
}

sub get_next_pos {
   my ($curr, $base, @seq_array) = @_;
   
   for (my $i = $curr; $i <= scalar(@seq_array); $i++ ) {
      if ( $base eq $seq_array[$i] ) {
         return ($i);
      }
   }
   print "Problem at $curr, base = $base\n";
   return("Problem at $curr");
}
