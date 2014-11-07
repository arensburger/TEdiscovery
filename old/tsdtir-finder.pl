#!/usr/bin/perl

# May 2014 finds positions of TIR-TSD combination in genome

use strict;
use File::Temp ();
use Getopt::Long;

my $MAXTIRLENGTH = 30;
my $MAXDISTANCE = 5000;
my $MINDISTANCE = 500;
my $genomefile;
my $outputname;
#####read and check the inputs
GetOptions(
	'g:s'   => \$genomefile,
	'o:s'	=> \$outputname
);
unless (($genomefile) and ($outputname)){
	die "usage: perl tsdtir-finder -g <fasta genome file, REQUIRED> -o <base name for bed outut file, <REQUIRED>\n";
}
my $bed1name = $outputname . "_elements.bed";
open (OUTPUT1, ">$bed1name") or die "cannot create file $bed1name\n";
my $bed2name = $outputname . "_tirs.bed";
open (OUTPUT2, ">$bed2name") or die "cannot create file $bed2name\n";

my %genome=genometohash($genomefile);
foreach my $segmentname (keys %genome) {
	my $segmentseq = File::Temp->new( UNLINK => 1, SUFFIX => '.fas' ); #sequence segment
	open (OUTPUT, ">$segmentseq") or die;
	print OUTPUT ">$segmentname\n";
	print OUTPUT "$genome{$segmentname}";
	close OUTPUT;
	my $einverted_out1 = File::Temp->new( UNLINK => 1, SUFFIX => '.blast' ); #sequence segment
	my $einverted_out2 = File::Temp->new( UNLINK => 1, SUFFIX => '.fas' ); #sequence segment
	`einverted -sequence $segmentseq -gap 12 -threshold 35 -match 3 -mismatch -4 -maxrepeat $MAXDISTANCE -outfile $einverted_out1 -outseq $einverted_out2 2>&1`;

	open (INPUT, $einverted_out2) or die;
	while (my $line = <INPUT>) {
		my $name1 = $line;
		my $seq1 = <INPUT>;
		my $name2 = <INPUT>;
		my $seq2 = <INPUT>;
		
		# check length
		if (length $seq1 < $MAXTIRLENGTH) { 
			my @data1 = split("_", $name1);
			my @data2 = split("_", $name2);
			my $b1 = $data1[-2]; # left bound or left TIR
			my $b2 = $data2[-1]; # right bound of right TIR

			if ($b1 > $b2) { #this happens if the lines are reversed
				my $temp = $b1;
				$b1 = $b2;
				$b2 = $temp;
				$temp = $seq1;
				$seq1 = $seq2;
				$seq2 = $temp;
			}

			next if (($b2 - $b1) < $MINDISTANCE);
			next if ($b2 - $b1) > $MAXDISTANCE;
	
			#check if has TA ends
			my $l1=substr($seq1, 0, 2);
			my $l2=substr($seq1, 1, 2);
			my $l3=substr($seq1, 2, 2);
			my $l4=substr($seq1, 3, 2);
			my $r1=substr($seq2, -3, 2);
			my $r2=substr($seq2, -4, 2);
			my $r3=substr($seq2, -5, 2);
			my $r4=substr($seq2, -6, 2);
			my $loffset=0;
			my $roffset=0;
			if ($l4 eq "ta") {
				$loffset=4;
			}
			elsif ($l3 eq "ta") {
				$loffset=3;
			}
			elsif ($l2 eq "ta") {
				$loffset=2;
			}
			elsif ($l1 eq "ta") {
				$loffset=1;
			}
			if ($r4 eq "ta") {
				$roffset=4;
			}
			elsif ($r3 eq "ta") {
				$roffset=3;
			}
			elsif ($r2 eq "ta") {
				$roffset=2;
			}
			elsif ($r1 eq "ta") {
				$roffset=1;
			}

			if ($loffset and $roffset) { # if offset means that TIRs are found
				$b1 = $b1 + $loffset;
				$b2 = $b2 - $roffset - 1;
				my $tsd1 = substr $seq1, $loffset + 1, (length $seq1) - $loffset -2;
				my $tsd2 = substr $seq2, 0, -($roffset+2);
				my $b1_2 = $b1 + (length $tsd1); # right bound of left TIR
				my $b2_2 = $b2 - (length $tsd2); 
				print OUTPUT1 "$segmentname\t$b1\t$b2\tTA-tir\n";
				print OUTPUT2 "$segmentname\t$b1\t$b1_2\t$segmentname\t$b2_2\t$b2\tTA-TIR\t0\t+\t+\ttir1=$tsd1;tir2=$tsd2\n";
			}
		}
	}	
}

#load a genome into a hash
sub genometohash {
	use strict;
	(my $filename) = @_;
	my %genome; #hash with the genome
	my $seq;
	my $title;
	open (INPUT, $filename) or die "cannot open input file $filename in sub genometohash\n";
	while (my $line = <INPUT>) {
		if (($line =~ />(\S+)/) && (length $seq > 1)) {
			if (exists $genome{$title}) {
				print STDERR "error in sub genometohash, two contigs have the name $title, ignoring one copy\n";
			}
			else {
				$genome{$title} = $seq;
			}
			#chomp $line;
			my @data = split(" ", $line);
			$title = $data[0];
			$title = substr $title, 1; # strip the > character from title
			$seq = "";
		}
		elsif ($line =~ />(\S+)/) { #will only be true for the first line
			#chomp $line;
			my @data = split(" ", $line);
			$title = $data[0];
			$title = substr $title, 1; # strip the > character from title
                        $seq = "";
		}
		else {
			$line =~ s/\s//g;
			$seq .= $line;
		}
	} 
	$genome{$title} = $seq;

	return (%genome);
}

