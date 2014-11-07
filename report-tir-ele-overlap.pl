#!/usr/bin/perl

use strict;
use Getopt::Long;

my $SIDE_SEQUENCE_LENGTH = 20; #number of bp to to get on each side of the element
#my $TSD_LENGTH = 2;
#my $TIR_LENGTH = 11;

my %config;
GetOptions(\%config,
	'rmbed=s',
	'tirbed=s',
	'genome=s',
	'out=s',
);

##> Check if no mandatory parameter is missing and set defaults
if (!exists $config{rmbed})         {printUsage();}
if (!exists $config{tirbed})         {printUsage();}
if (!exists $config{genome})         {printUsage();}
if (!exists $config{out}) {$config{out} = "out";}

## Main program ##
my %species; # holds the name of the element as key and an array with each instance detail as value 
### find overlaps between the databases using bedtools ###
#my $overlapfile = `bedtools intersect -a $config{tirbed} -b $config{rmbed} -wo`;
my @lines = split /\n/, `bedtools intersect -a $config{tirbed} -b $config{rmbed}  -wo`; 
unless (@lines) { # checking that something got put into the array
	die "ERROR, bedtools intersection did not work, check that bedtools is installed\n";
}


### read through the file ###
my %seen; #to prevent duplications
foreach my $line (@lines) {
	my @data = split " ", $line;
	my $scaffold = $data[0];
	my $tirb1 = $data[1];
	my $tirb2 = $data[2];
	my $eleb1 = $data[5];
	my $eleb2 = $data[6];
	my $type = $data[7];
	my $ori = $data[9];
	my $overlap = $data[10];
	#determine the TSD and TIR length for current element
	my $tsd_length;
	my $tir_length;
	if($data[3] =~ /(\S+)_(\S+)/) {
		my $tsd_data = $1;
		my $tir_data = $2;

		if ($tsd_data eq "TA") {
			$tsd_length = 2;
		}
		elsif ($tsd_data =~ /(\d+)TSD/) {
			$tsd_length = $1;
		}
		else {
			die "don't know TSD $tsd_data from line\n$line\n";
		}

		if ($tir_data =~ /(\d+)TIR/) {
			$tir_length = $1;
		}
		else {
			die "don't know TIR $tir_length from line\n$line\n";
		}
	
	}	

	if (($eleb2 - $eleb1) == $overlap) { #true if the sequences overlap fully
		my @data2 = split ("#", $type);
#		print "$data2[0], $data2[1], $data2[2], $type\n";
#		print "$scaffold, $tirb1, $tirb2, $eleb1, $eleb2, $type\n";
		# add line into %species
		my $seen_summary =  "$scaffold\t$tirb1\t$tirb2\t$ori\t$data[7]";
		unless (exists $seen{$seen_summary}) {
			push @{ $species{$data2[0]} }, "$scaffold\t$tirb1\t$tirb2\t$eleb1\t$eleb2\t$ori\t$data2[1]\t$data2[2]\t$tsd_length\t$tir_length\n";
			$seen{$seen_summary} = 0;
		}
	}
}

### report data on each element ###
my %genome = genometohash($config{genome});
for my $element ( keys %species ) {
	my @element_data; # holds all the elements for this species
#print "$element:\n";
	for my $i ( 0 .. $#{ $species{$element} } ) {
		my @data = split("\t", $species{$element}[$i]);

		#get the tir and tsd lengths
		my $tsd_length = $data[8];
		my $tir_length = $data[9];

		#get the element sequence
		my $element_sequence;
		if ($data[5] eq "+") {
			$element_sequence = substr($genome{$data[$data[0]]}, ($data[1] - $SIDE_SEQUENCE_LENGTH), ($data[2] - $data[1] + $SIDE_SEQUENCE_LENGTH + $SIDE_SEQUENCE_LENGTH)); 
		}
		elsif ($data[5] eq "-") {
			$element_sequence = rc(substr($genome{$data[$data[0]]}, ($data[1] - $SIDE_SEQUENCE_LENGTH), ($data[2] - $data[1] + $SIDE_SEQUENCE_LENGTH + $SIDE_SEQUENCE_LENGTH))); 
		}
		else {
			die "unexpected orientation $data[5]\n";
		}
	
		my $left_side_sequence = substr($element_sequence, 0, $SIDE_SEQUENCE_LENGTH);
		my $left_TSD = substr($element_sequence, $SIDE_SEQUENCE_LENGTH, $tsd_length);
		my $left_TIR = substr($element_sequence, ($SIDE_SEQUENCE_LENGTH + $tsd_length), $tir_length);
		my $right_side_sequence = substr($element_sequence, -$SIDE_SEQUENCE_LENGTH, $SIDE_SEQUENCE_LENGTH);
		my $right_TSD = substr($element_sequence, -(length($right_side_sequence) + $tsd_length), $tsd_length);
		my $right_TIR = substr($element_sequence, -($SIDE_SEQUENCE_LENGTH + $tsd_length + $tir_length), $tir_length); 
		my $te_sequence = substr($element_sequence, ($SIDE_SEQUENCE_LENGTH + $tsd_length), (length $element_sequence) - $SIDE_SEQUENCE_LENGTH - $SIDE_SEQUENCE_LENGTH - $tsd_length - $tsd_length); 

#print "$element_sequence\n";
#print "$left_side_sequence\t$left_TIR\t$right_TIR\t$right_side_sequence\t$data[0]\t$data[1]\t$data[2]\n";
#print "$left_TSD, $right_TSD\n";
#print "$te_sequence\n";
		# store this element, the test prevents adding elements that could not be properly determined
		unless ((length $left_side_sequence == 0) or (length $left_TSD == 0) or (length $left_TIR == 0)) {
			my $location = $data[0] . ":" . $data[1] . "-" . $data[2];
			chomp $data[7];
			my $id = $data[6] . "/" . $data[7];
			push (@element_data, "$left_side_sequence\t$left_TIR\t$right_TIR\t$right_side_sequence\t$location\t$id\t$te_sequence\t$left_TSD\t$right_TSD");
		}
	}

	#compare all the elements for this species to each other
	my %insert_data; #holds all the information about insertion data for each element
	my %sequence_data; #holds the sequence for each element

	for (my $i=0; $i< ((scalar @element_data) - 1); $i++) {
		my @data = split("\t", $element_data[$i]);
		my $location1 = $data[4];
		my $tirs1 = $data[1] . $data[2];
		my $side_sequences1 = $data[0] . $data[3];
		my $id1 = $data[5];
		my $tsd1 = $data[7];
		my $te_sequence1 = ">" . $data[4] . "#" . $id1 . "#" . $element . "\n" . $data[6];
		for (my $j=$i+1; $j<scalar @element_data; $j++) {
			my @data2 = split("\t", $element_data[$j]);
			my $location2 = $data2[4];
			my $tirs2 = $data2[1] . $data2[2];
			my $side_sequences2 = $data2[0] . $data2[3];
			my $id2 = $data2[5];
			my $tsd2 = $data2[8];
			my $te_sequence2 = ">" . $data2[4] . "#" . $id2 . "#" . $element . "\n" . $data2[6];

			# testing for uniqueness
			my $TIR_IDENTITY_THRESHOLD = 0.1; # maximum p distance to call two tirs identical
			my $SIDE_IDENTITY_THRESHOLD = 0.7; # minimum distance between the side sequences to call them different 
			if ($location1 ne $location2) {
				my $passed_tir_test = 0; # boolean, 0 if tirs are identical enough 1 if not
				my $passed_side_test = 0; # boolean, 0 if side sequence are the same 1 if they are different
				if (pdistance ($tirs1, $tirs2) < $TIR_IDENTITY_THRESHOLD) {
					$passed_tir_test = 1;
				}
				if (pdistance ($side_sequences1, $side_sequences2) > $SIDE_IDENTITY_THRESHOLD) {
					$passed_side_test = 1;
				}

				if ($passed_tir_test and $passed_side_test) {
					my $insertion_data1 = "$data[0]\t$data[7]\t$data[1]\t$data[2]\t$data[8]\t$data[3]\t$location1";
					my $insertion_data2 = "$data2[0]\t$data2[7]\t$data2[1]\t$data2[2]\t$data2[8]\t$data2[3]\t$location2";
					$insert_data{$insertion_data1} = "$data[1]\t$data[2]";
					$insert_data{$insertion_data2} = "$data2[1]\t$data2[2]";
					$sequence_data{$te_sequence1} = 0;
					$sequence_data{$te_sequence2} = 0;
				}
			}			
		}
	}

	# print if any results found
	if (%insert_data) {
		print "********\n";
		foreach my $key (sort {$insert_data{$a} cmp $insert_data{$b} } keys %insert_data) {
			print "$key\n";
		}
		foreach my $key (keys %sequence_data) {
			print "$key\n";
		}
		print "\n";
	}
}

sub printUsage{

print STDOUT "DESCRIPTION: This program takes bed formated files from RepeatMasker and from genome TIR analysis and returns report\n";
print STDOUT "USAGE : combine-tir -r \"RepeatMasker .bed file\" -t \"TIR analysis .bed file\" -g \"Genome fasta file\" -o \"output file\"
Options : 
    -r | rmbed		RepeatMasker output from RepeatExplorer file, bed formated (Mandatory)
    -t | tirbed		Whole genome TIR analysis, bed formated (Mandatory)
    -g | genome		Fasta formated genome (Mandatory)
    -o | out   		Name of ouput file (default \"out\")\n";
    exit;
}

#load a genome into a hash
sub genometohash {
	use strict;
	(my $filename) = @_;
	my %genome; #hash with the genome
	my $seq="";
	my $title;
	open (INPUT, $filename) or die "cannot open input file $filename in sub genometohash\n";
	while (my $line = <INPUT>) {
		if (($line =~ />(\S+)/) && (length $seq > 1)) {
			if (exists $genome{$title}) {
				print STDERR "error in sub genometohash, two contigs have the name $title, ignoring one copy\n";
#				exit;
			}
			else {
				$genome{$title} = $seq;
			}
			#chomp $line;
			my @data = split(" ", $line);
			$title = $data[0];
			$title = substr $title, 1;
			$seq = "";
		}
		elsif ($line =~ />(\S+)/) { #will only be true for the first line
			#chomp $line;
			my @data = split(" ", $line);
			$title = $data[0];
		$title = substr $title, 1;
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
sub rc {
    my ($sequence) = @_;
    $sequence = reverse $sequence;
    $sequence =~ tr/ACGTRYMKSWacgtrymksw/TGCAYRKMWStgcayrkmws/;
    return ($sequence);
}

# takes two sequences and returns distance value
sub pdistance {
	(my $seq1, my $seq2) = @_;
	my $differences;
	if (length $seq1 == length $seq2) {
		for (my $i=0; $i<(length $seq1); $i++) {
			if ((uc(substr $seq1, $i, 1)) ne (uc(substr $seq2, $i, 1))) {
				$differences++;
			}
		}
	}
	else {
		die "different sequence lengths\n";
	}

	my $pdistance = $differences/(length $seq1);
	return($pdistance);
}
exit;

