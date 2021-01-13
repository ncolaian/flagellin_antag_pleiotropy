#! /usr/bin/env perl

# This will create a matrix that can be used for MIT

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Data::Dumper;
use Readonly;
use YAML qw( DumpFile LoadFile );

#My Variables
my $help = 0;
my $man = 0;
my $aa;
my $msa;
my $out;
my $matrix_builder = "/pine/scr/n/c/ncolaian/mamp_pa/bin/mutual_info_2/MIT_matrix_builder_part.pl";
my $mit_calculator = "/pine/scr/n/c/ncolaian/mamp_pa/bin/mutual_info_2/MIT_calculator.pl";
my $mit_type = 2; #2 was shown to work best on simulated data. Use this

#Read in the variables from the command line
GetOptions( 'man'   =>  \$man,
            'help'  =>  \$help,
            'aa'    =>  \$aa,
            'msa|m=s'   =>  \$msa,
            'out|o=s'   =>  \$out,
            'mit_type|t:i'  =>  \$mit_type,
            ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manal and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

Readonly::Array my @AMINO_ACIDS => qw/
A R N D C Q E G H I L K M F P S T W Y V
/;

Readonly::Array my @N_ACIDS => qw/ A T G C /;

#main();

if ( $mit_type == 2 ) {
    main2();
}
if ( $mit_type == 3 ) {
    main3();
}

## Subroutines ##
sub main {
    my $seq_aref = get_sequences();
    my $pos_href = get_hash_pos($seq_aref);
    create_matrix($pos_href);
    combine_matrix_parts($pos_href);
    my $pos_sum_href = pos_sums_per_pos($pos_href);
    perform_MIT($pos_sum_href);
}

sub main2{
    my $seq_aref = get_sequences();
    my $pos_href = get_hash_pos($seq_aref);
    create_matrix($pos_href);
    combine_matrix_parts($pos_href);
    perform_MIT_type2();
}

sub main3{
    my $seq_aref = get_sequences();
    my $pos_href = get_hash_pos($seq_aref);
    create_matrix($pos_href);
    combine_matrix_parts($pos_href);
    perform_MIT_type3();
}


sub get_sequences {
    open my $IN, "<", $msa;
    my $thr_awy = <$IN>;
    
    my @seq;
    my $mm_seq = "";
    while ( <$IN> ) {
        chomp $_;
        if ( $_ =~ /^>/ ) {
            #handle first line
            if ( $mm_seq ne "" ) {
                push @seq, $mm_seq;
                $mm_seq = "";
                next;
            }
            next;
        }
        $mm_seq .= $_;
    }
    #handle last line
    push @seq, $mm_seq;
    close($IN);
    return( \@seq );
}

sub get_hash_pos {
    my ( $seq_aref ) = @_;
    
    my %pos_hash;
    foreach my $seq ( @$seq_aref ) {
        my @bases = split //, $seq;
        
        my $pos = 1;
        foreach my $b ( @bases ) {
            if ( $pos_hash{$pos} ) {
                push @{$pos_hash{$pos}}, $b;
            }
            else {
                $pos_hash{$pos} = [$b];
            }
            $pos++;
        }
    }
    
    return( \%pos_hash );
}

sub create_matrix {
    my ( $pos_href ) = @_;
    
    #create a file containing the hash
    my @jobs;
    DumpFile("$out/build.hash", $pos_href);
    
    my $num_pos = keys %$pos_href;
    my $interval = int($num_pos/20);
#    for ( my $i = 1; $i <= $num_pos; $i = $i+$interval ) {
 #       my $end = $i + $interval - 1;
  #      #handles potential over flow
   #     if ( $end > $num_pos ) {
    #        $end = $num_pos;
     #   }
      #  my $cmd;
       # if ( $aa ) {
        #    $cmd = "perl $matrix_builder $i $end a $out";
        #}
        #else {
        #    $cmd = "perl $matrix_builder $i $end n $out"
        #}
        #push @jobs,$cmd;
    #}
    
    #submit_and_stall_slurm(\@jobs, "matrix_parts", $out, "4000", "1:00:00");
    `perl $matrix_builder 1 $num_pos a $out`;
    return(1);
}

sub combine_matrix_parts {
    my ( $pos_href ) = @_;
    open my $OUT, ">", "$out/mit_mat_full.txt";
    
    #print header
    print $OUT "i\tj";
    if ( $aa ) {
        foreach my $aa ( @AMINO_ACIDS ) {
            foreach my $aa2 ( @AMINO_ACIDS ) {
                my $key = "$aa-$aa2";
                print $OUT "\t$key";
            }
        }
    }
    else {
        foreach my $aa ( @N_ACIDS ) {
            foreach my $aa2 ( @N_ACIDS ) {
                my $key = "$aa-$aa2";
                print $OUT "\t$key";
            }
        }
    }
    print $OUT "\n";
    close $OUT;
    
    #combine the parts
    my $num_pos = keys %$pos_href;
    my $interval = int($num_pos/20);
    for ( my $i = 1; $i <= $num_pos; $i = $i+$interval ) {
        `cat $out/mat_pt$i.txt >> $out/mit_mat_full.txt`;
        `rm $out/mat_pt$i.txt`;
    }
    return(1);
}

sub pos_sums_per_pos {
    my ($href) = @_;
    
    my %sums;
    foreach my $pos ( keys %$href ) {
        my %hash;
        my $bases = $href->{$pos};
        
        my $total = 0;
        foreach my $b ( @$bases ) {
            next if ( $b eq "-");
            $total++;
            if ( $hash{$b} ) {
                $hash{$b} += 1;
            }
            else {
                $hash{$b} = 1;
            }
        }
        $hash{'total'} = $total;
        $sums{$pos} = \%hash;
    }
    
    return(\%sums);
}

sub perform_MIT_para {
    my ( $sums_href ) = @_;
    DumpFile("$out/sums.hash", $sums_href);
    
    my @jobs;
    my $mit_calcs = `wc -l $out/mit_mat_full.txt`;
    $mit_calcs =~ int((split /\s/, $mit_calcs)[0]);
    
    my $interval = int($mit_calcs/20);
    for ( my $i = 1; $i <= $mit_calcs; $i = $i+$interval ) {
        my $end = $i + $interval - 1;
        if ( $end > $mit_calcs ) {
            $end = $mit_calcs;
        }

        my $cmd = "perl $mit_calculator $i $end $out";

        push @jobs,$cmd;
    }
    
    submit_and_stall_slurm(\@jobs, "matrix_calc", $out, "400", "1:00:00");
    
    #Combine all the parts
    open my $OUT, "<", "$out/mit.txt";
    print $OUT "Pos1\tPos2\tEntropy\n";
    close $OUT;
    for ( my $i = 1; $i <= $mit_calcs; $i = $i+$interval ) {
        `cat $out/calc_pt$i.txt >> $out/mit.txt`;
        `rm $out/calc_pt$i.txt`;
    }
    return(1);
}

sub perform_MIT {
    my ( $sums_href ) = @_;
    open my $MATRIX, "<", "$out/mit_mat_full.txt";
    open my $OUT, ">", "$out/mit.txt";
    print $OUT "Pos1\tPos2\tEntropy\n";
    
    my $header = <$MATRIX>;
    chomp $header;
    my @header_split = split /\t/, $header;
    
    my %pos;
    for ( my $i = 2; $i < scalar(@header_split); $i++ ) {
        my @mm_array = (split /-/, $header_split[$i]);
        $pos{$i} = \@mm_array;
    }
    
    my %already_done;
    while ( <$MATRIX> ) {
        chomp $_;
        print "$_\n";
        my @split_line = split /\t/, $_;
        #if the comparison is to itself skip MIT
        next if ( $split_line[0] == $split_line[1] );
        my $pos1 = shift @split_line;
        my $pos2 = shift @split_line;
        
        #only perform MIT once between base pair positions
        if ( $already_done{$pos2} ) {
            next;
        }
        else {
            if ( !($already_done{$pos1}) ) {
                $already_done{$pos1} = 1;
            }
        }
        
        #get pos hrefs
        my $pos1_sum_href = $sums_href->{$pos1};
        my $pos2_sum_href = $sums_href->{$pos2};
        
        #now need to perform MIT
        my $total_comb = array_sum(\@split_line);
        if ( $total_comb == 0 ) {
            next;
        }
        
        my $total_mit = 0.000000;
        for ( my $i = 0; $i < scalar(@split_line); $i++ ) {
            next if ( $split_line[$i] == 0 );
            
            my @bases = @{$pos{$i+2}}; #go through the paired amino acids
            my $prob_1 = 0;
            my $prob_2 = 0;
            eval {
                $prob_1 = $pos1_sum_href->{$bases[0]}/$pos1_sum_href->{'total'};
                $prob_2 = $pos2_sum_href->{$bases[1]}/$pos2_sum_href->{'total'};
            } or do {
                next;
            };
            
            my $comb_prob = $split_line[$i]/$total_comb;
            if ( $comb_prob == 0 || $prob_1 == 0 || $prob_2 == 0) {
                next;
            }
            
            
            $total_mit += $comb_prob * log($comb_prob/($prob_1*$prob_2));
        }
        print $OUT "$pos1\t$pos2\t$total_mit\n";
    }
    close $MATRIX;
    close $OUT;
}

sub perform_MIT_type2 {
    open my $MATRIX, "<", "$out/mit_mat_full.txt";
    open my $OUT, ">", "$out/mit.txt";
    print $OUT "Pos1\tPos2\tEntropy\n";
    
    my $header = <$MATRIX>;
    chomp $header;
    my @header_split = split /\t/, $header;
    
    my %pos;
    for ( my $i = 2; $i < scalar(@header_split); $i++ ) {
        my @mm_array = (split /-/, $header_split[$i]);
        $pos{$i} = \@mm_array;
    }
    
    my %already_done;
    while ( <$MATRIX> ) {
        chomp $_;
        my @split_line = split /\t/, $_;
        #if the comparison is to itself skip MIT
        next if ( $split_line[0] == $split_line[1] );
        my $pos1 = shift @split_line;
        my $pos2 = shift @split_line;
        
        #only perform MIT once between base pair positions
        if ( $already_done{$pos2} ) {
            next;
        }
        else {
            if ( !($already_done{$pos1}) ) {
                $already_done{$pos1} = 1;
            }
        }
        
        #start 2 counters for p(ai) and p(bj)
        my %pai;
        my %pbj;
        if ( $aa ) {
            foreach my $base ( @AMINO_ACIDS ) {
                $pai{$base} = 0;
                $pbj{$base} = 0;
            }
        }
        else {
            foreach my $base ( @N_ACIDS ) {
                $pai{$base} = 0;
                $pbj{$base} = 0;
            }
        }
        
        #now need to perform MIT
        my $total_comb = array_sum(\@split_line);
        if ( $total_comb == 0 ) {
            next;
        }
        
        my $total_mit = 0.000000;
        #first go through and get the per base probability
        for ( my $i = 0; $i < scalar(@split_line); $i++ ) {
            next if ( $split_line[$i] == 0 );
            
            my @bases = @{$pos{$i+2}}; #go through the paired amino acids
            
            my $comb_prob = $split_line[$i]/$total_comb;
            $pai{$bases[0]} += $comb_prob;
            $pbj{$bases[1]} += $comb_prob;
        }
        for ( my $i = 0; $i < scalar(@split_line); $i++ ) {
            next if ( $split_line[$i] == 0 );
            
            my @bases = @{$pos{$i+2}}; #go through the paired amino acids
            
            my $comb_prob = $split_line[$i]/$total_comb;
            #perform the math
            $total_mit += $comb_prob*log( $comb_prob/($pai{$bases[0]} * $pbj{$bases[1]}) );
        }
        
        print $OUT "$pos1\t$pos2\t$total_mit\n";
    }
    close $MATRIX;
    close $OUT;
}

#This is from martin et al. 2005
sub perform_MIT_type3 {
    open my $MATRIX, "<", "$out/mit_mat_full.txt";
    open my $OUT, ">", "$out/mit.txt";
    print $OUT "Pos1\tPos2\tEntropy\n";
    
    my $header = <$MATRIX>;
    chomp $header;
    my @header_split = split /\t/, $header;
    
    my %pos;
    for ( my $i = 2; $i < scalar(@header_split); $i++ ) {
        my @mm_array = (split /-/, $header_split[$i]);
        $pos{$i} = \@mm_array;
    }
    
    my %already_done;
    while ( <$MATRIX> ) {
        chomp $_;
        my @split_line = split /\t/, $_;
        #if the comparison is to itself skip MIT
        next if ( $split_line[0] == $split_line[1] );
        my $pos1 = shift @split_line;
        my $pos2 = shift @split_line;
        
        #only perform MIT once between base pair positions
        if ( $already_done{$pos2} ) {
            next;
        }
        else {
            if ( !($already_done{$pos1}) ) {
                $already_done{$pos1} = 1;
            }
        }
        
        #start 2 counters for p(ai) and p(bj)
        my %pai;
        my %pbj;
        if ( $aa ) {
            foreach my $base ( @AMINO_ACIDS ) {
                $pai{$base} = 0;
                $pbj{$base} = 0;
            }
        }
        else {
            foreach my $base ( @N_ACIDS ) {
                $pai{$base} = 0;
                $pbj{$base} = 0;
            }
        }
        
        #now need to perform MIT
        my $total_comb = array_sum(\@split_line);
        if ( $total_comb == 0 ) {
            next;
        }
        
        my $total_mit = 0.000000;
        
        my %combined_prob;
        #first go through and get the per base entropy
        for ( my $i = 0; $i < scalar(@split_line); $i++ ) {
            next if ( $split_line[$i] == 0 );
            
            my @bases = @{$pos{$i+2}}; #go through the paired amino acids
            
            my $comb_prob = $split_line[$i]/$total_comb;
            $pai{$bases[0]} += $comb_prob;
            $pbj{$bases[1]} += $comb_prob;
            $combined_prob{$split_line[$i]} = $comb_prob;
        }
        ############################
        ## H(a)+H(b)-H(a,b) -> MI ##
        ############################
        ###
        ### I have also changed the scaling log factor due to Oliviera et al. 2007
        ###
        #calculate H(a) and H(b)
        my $a_entropy = 0.00000;
        my $b_entropy = 0.00000;
        foreach my $base ( keys %pai ) {
            if ( $pai{$base} != 0 ){
                my $value = $pai{$base} * ( log($pai{$base}));#/log($total_comb) );
                $a_entropy += $pai{$base} * ( log($pai{$base}));#/log($total_comb) );
            }
            if ( $pbj{$base} != 0 ) {
                $b_entropy += $pbj{$base} * ( log($pbj{$base}));#/log($total_comb) );
            }
            #we devide the natural log by the log of the number we want the base log to be
            #this scales the entropy to 1
        }
        #calculate H(a,b)
        my $combined_entropy = 0.00000;
        foreach my $part ( keys %combined_prob ) {
            if ( $combined_prob{$part} != 0 ) {
                #to get a log base to the correct base you do the log(number)/log(base wanted)
                $combined_entropy += $combined_prob{$part} * ( log($combined_prob{$part}));#/log($total_comb) );
            }
        }
        
        #calculate mit
        $total_mit = -($a_entropy) + -($b_entropy) - -($combined_entropy);
        
        #scale mit by the combined entropy
        if ( $combined_entropy != 0 ) {
            $total_mit = $total_mit/-($combined_entropy);
        }
        
        print $OUT "$pos1\t$pos2\t$total_mit\n";
    }
    close $MATRIX;
    close $OUT;
}

sub array_sum {
    my ($aref) = @_;
    my $sum = 0;
    foreach my $part ( @$aref) {
        $sum += $part;
    }
    return $sum;
}

sub submit_and_stall_slurm {
    my ( $job_aref, $job_name, $out_directory, $mem, $time ) = @_;
    $logger->info("Submitting and Stalling until all $job_name jobs finish\n");
    
    foreach my $cmd ( @$job_aref ) {
        my $sbatch_cmd = "sbatch -t $time -o $out_directory/$job_name.out --open-mode=append --mem=$mem -J $job_name --wrap=\"$cmd\"";
        system($sbatch_cmd);
    }
    
    #Create Active job file
    `echo start > $out_directory/ACTIVE_JOBS`;
    
    my $minutes = 0;
    while ( -s "$out_directory/ACTIVE_JOBS" ) {
        sleep 60;
        `squeue --name=$job_name -h > $out_directory/ACTIVE_JOBS`;
        $minutes++;
    }
    
    `rm $out_directory/ACTIVE_JOBS`;
    $logger->info("All $job_name jobs have finished in ~$minutes minutes");
    
    return 1;
}

__END__
=head1 TITLE



=head1 VERSION



=head1 INCLUDED MODULES

Getopt::Long;
Pod::Usage;
Carp;
Readonly;
Path::Class;
Data::Dumper;
Log::Log4perl qw(:easy);
Log::Log4perl::CommandLine qw(:all);

=head1 INHERIT

=head1 SYNOPSIS

=head1 PARAMETERS

=head1 CONFIGURATION AND ENVIRONMENT

    Need to load muscle and BLAST

=head1 DEPENDENCIES

    muscle
    BLAST
    Master_aln
    
=head1 INCOMPATIBILITIES

    None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests	
	
=head1 AUTHOR

Nicholas Colaianni
contact via C<< <ncolaian@live.unc.edu> >>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2017, Nicholas Colaianni
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut
