#! /usr/bin/env perl

# This script will create a random file containing a number of sequence and combinations used to test MIT

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

#My Variables
my $help = 0;
my $man = 0;
my $out;
my $seqs;
my $null = 0;
my $real_null;

#Read in the variables from the command line
GetOptions( 'man'   =>  \$man,
            'help'  =>  \$help,
            'out|o=s'   =>  \$out,
            'seqs|s=i'  =>  \$seqs,
            'null|n:i'  =>  \$null,
            'real_n|r:i'    =>  \$real_null,
            ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manal and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

#array for big random
Readonly::Array my @AMINO_ACIDS => qw/
A R N D C Q E G H I L K M F P S T W Y V
/;

## Main ##
open my $OUT, ">", $out;
#sequences will be 20 amino acids long
for (my $i = 0; $i < $seqs; $i++ ){
    my @array;
    
    if ( $real_null ) {
        $array[$real_null-1] = "A";
        @array = ("A") x @array;
        my $changes = 1+int(rand(2));
        my @changes_ar;
        $changes_ar[$changes-1] = 0;
        foreach my $c ( @changes_ar ) {
            $c = int(rand($real_null));
            $array[$c] = $AMINO_ACIDS[rand @AMINO_ACIDS];
        }
        print $OUT ">Seq$i\n" . join("", @array) . "\n";
    }
    
    elsif ( $null ) {
        $array[$null-1] = "A";
        for ( my $i = 0; $i < $null; $i++ ) {
            $array[$i] = $AMINO_ACIDS[rand @AMINO_ACIDS];
        }
        print $OUT ">Seq$i\n" . join("", @array) . "\n";
    }
    else {
        $array[19] = "A";
        @array = ("A") x @array;
        
        #amino acids at position 2 and 19 will be linked to change the same time to the same random amino acid
        my $rand = $AMINO_ACIDS[rand @AMINO_ACIDS];
        $array[1] = $rand;
        $array[18] = $rand;
        
        #amino acids 4 and 17 will both change every time but to random amino acids
        $array[3] = $AMINO_ACIDS[rand @AMINO_ACIDS];
        $array[16] = $AMINO_ACIDS[rand @AMINO_ACIDS];
        
        #amino acids 6 and 15 will change randomly at about 50% of the time to the same amino acid
        if ( rand() > .5 ) {
            $rand = $AMINO_ACIDS[rand @AMINO_ACIDS];
            $array[5] = $rand;
            $array[14] = $rand;
        }
        
        #amino acids 8 and 13 will change randomly at about 50% of the time to random amino acids
        if ( rand() > .5 ) {
            $array[7] = $AMINO_ACIDS[rand @AMINO_ACIDS];
            $array[12] = $AMINO_ACIDS[rand @AMINO_ACIDS];
        }
        
        #amino acids 10 and 11 will change 50% and to the same amino acid every time
        if ( rand() > .5 ) {
            $array[9] = "N";
            $array[10] = "N";
        }
        print $OUT ">Seq$i\n" . join("", @array) . "\n";
    }
}

close($OUT);



## Subroutines ##

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