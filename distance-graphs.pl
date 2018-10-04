#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Temp ();

my $RMfile; 
my $genome;
my $outputfile;
#####read and check the inputs
GetOptions(
	'r:s'   => \$RMfile,
	'g:s'	=> \$genome,
	'o:s'	=> \$outputfile
);
unless ($RMfile and $genome and $outputfile) {
	die "usage: perl distance-graph.pl -r <RepeatMasker file, REQUIRED> -g <genome file, REQUIRED> -o <output file, REQUIRED>";
}

# read the RM file and create a bed file of relevant sequences
my $rmbedfile = File::Temp->new( UNLINK => 1, SUFFIX => '.bed' );
open (OUTPUT, ">$rmbedfile") or die "cannot open $rmbedfile for writting";
open (INPUT, $RMfile) or die "cannot open file $ARGV[0]\n";
<INPUT>;
<INPUT>;
<INPUT>;
while (my $line = <INPUT>) {
	my @data = split (" ", $line);
	my $class; # TE class
	if ($data[10] =~ /(\S+)\/(\S+)/) { # only read if the class and subfamily are specfied
		$class = $1;
		

	}
}

my $seq1 = "AACGTGGCCACAT";
my $seq2 = "AAGGTCGCCACAC";
my $distance = k2p($seq1, $seq2);
print "$distance\n";

sub k2p {
	my ($seq1, $seq2) = @_;
	unless (length $seq1 == length $seq2) {
		return("NAL");
	}

	my @s1 = split("", uc($seq1));
	my @s2 = split("", uc($seq2));

	my $P; #number of transition positions
	my $Q; #number of transversion positions
	my $len; #length of compared sequence
	for (my $i=0; $i<scalar @s1; $i++) {
		if (("ATGC" =~ /$s1[$i]/) and ("ATGC" =~ /$s2[$i]/)) {
			$len++;
			unless ($s1[$i] eq $s2[$i]) {
				if ((($s1[$i] eq "A") and ($s2[$i] eq "G")) or (($s1[$i] eq "G") and ($s2[$i] eq "A")) or
				   (($s1[$i] eq "C") and ($s2[$i] eq "T")) or (($s1[$i] eq "T") and ($s2[$i] eq "C"))) {
					$P++;
				}
				else {
					$Q++;
				}
				
			} 
		}
	}
	if ($len) {
		my $Pprop = $P/$len;
		my $Qprop = $Q/$len;
		my $distance = (-0.5) * (log( (1-(2*$Pprop)-$Qprop) * sqrt(1-(2*$Qprop) )));
		return($distance);
	}
	else {
		return("NA0");
	}
}
