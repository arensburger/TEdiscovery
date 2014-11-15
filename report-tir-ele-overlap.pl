#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Temp();

my $SIDE_SEQUENCE_LENGTH = 20; #number of bp to to get on each side of the element

my %config;
GetOptions(\%config,
	'rmbed=s',
	'tirbed=s',
	'genome=s',
	'dbte=s',
	'out=s',
);

##> Check if no mandatory parameter is missing, set defaults, and output file
if (!exists $config{rmbed})         {printUsage();}
if (!exists $config{tirbed})         {printUsage();}
if (!exists $config{genome})         {printUsage();}
if (!exists $config{dbte})         {$config{dbte} = "./TEdiscovery/repbase_classII.fa";}
if (!exists $config{out}) {$config{out} = "out";}
open (OUTPUT, ">$config{out}") or die "cannot open output file $config{out}\n";
print OUTPUT "side LEFT\tTSD LEFT\tTIR LEFT\tTIR RIGHT\tTSD RIGHT\tside RIGHT\tBLAST name\tBLAST evalue\tBLAST overlap\tlength\tTIR length\tORF length\tsequence NAME\tsequence string\n";

## Main program ##

my %species; # holds the name of the element as key and an array with each instance detail as value 

### find overlaps between the databases using bedtools ###
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

	# determine the TSD and TIR length for current element
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

	# if the tir and predicted te sequences overlap (i.e. the te is inside the tirs) then record this event
	if (($eleb2 - $eleb1) == $overlap) { #true if the sequences overlap fully
		my @data2 = split ("#", $type);

		# add line into %species
		my $seen_summary =  "$scaffold\t$tirb1\t$tirb2\t$ori\t$data[7]";
		unless (exists $seen{$seen_summary}) {
			push @{ $species{$data2[0]} }, "$scaffold\t$tirb1\t$tirb2\t$eleb1\t$eleb2\t$ori\t$data2[1]\t$data2[2]\t$tsd_length\t$tir_length\n";
			$seen{$seen_summary} = 0;
		}
	}
}

### report data on each element ###

my %genome = uc_genometohash($config{genome}); # get genome and convert all letters to uppercase

for my $element (keys %species ) { # scroll throught elements
	my @element_data; # holds all the elements for this species of TE
	for my $i ( 0 .. $#{ $species{$element} } ) {
		my @data = split("\t", $species{$element}[$i]);

		## get the tir and tsd lengths
		my $tsd_length = $data[8];
		my $tir_length = $data[9];

		## get the element sequence
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
	
		## collect the rest of the information about this element 
		my $left_side_sequence = substr($element_sequence, 0, $SIDE_SEQUENCE_LENGTH);
		my $left_TSD = substr($element_sequence, $SIDE_SEQUENCE_LENGTH, $tsd_length);
		my $left_TIR = substr($element_sequence, ($SIDE_SEQUENCE_LENGTH + $tsd_length), $tir_length);
		my $right_side_sequence = substr($element_sequence, -$SIDE_SEQUENCE_LENGTH, $SIDE_SEQUENCE_LENGTH);
		my $right_TSD = substr($element_sequence, -(length($right_side_sequence) + $tsd_length), $tsd_length);
		my $right_TIR = substr($element_sequence, -($SIDE_SEQUENCE_LENGTH + $tsd_length + $tir_length), $tir_length); 
		my $te_sequence = substr($element_sequence, ($SIDE_SEQUENCE_LENGTH + $tsd_length), (length $element_sequence) - $SIDE_SEQUENCE_LENGTH - $SIDE_SEQUENCE_LENGTH - $tsd_length - $tsd_length); 

		## store the information into the @element_data array

		# this test prevents adding element data that could not be properly determined for one reason or another
		unless ((length $left_side_sequence == 0) or (length $left_TSD == 0) or (length $left_TIR == 0)) {
			my $location = $data[0] . ":" . $data[1] . "-" . $data[2];
			chomp $data[7];
			my $id = $data[6] . "/" . $data[7];
			push (@element_data, "$left_side_sequence\t$left_TIR\t$right_TIR\t$right_side_sequence\t$location\t$id\t$te_sequence\t$left_TSD\t$right_TSD");
		}
	}

	## compare all the elements for this TE species to each other

	my %insert_data; #holds all the information about insertion data for each element
	my %sequence_data; #holds the sequence for each element

	for (my $i=0; $i< ((scalar @element_data) - 1); $i++) { ## first loop 

		my @data = split("\t", $element_data[$i]);
		my $location1 = $data[4]; # genomic location of first member of pair
		my $tirs1 = $data[1] . $data[2]; # combined TIR sequences of first member of pair
		my $side_sequences1 = $data[0] . $data[3]; # combined side sequences of first member of pair
		my $id1 = $data[5]; # unique id of first member of pair
		my $tsd1 = $data[7]; # TSD sequence of first member of pair
		my $te_sequence1 = ">" . $data[4] . "#" . $id1 . "#" . $element . "\n" . $data[6]; # fasta formated genomic sequence of first member of pair

		for (my $j=$i+1; $j<scalar @element_data; $j++) { ## second loop

			my @data2 = split("\t", $element_data[$j]); 
			my $location2 = $data2[4]; # genomic location of second member of pair
			my $tirs2 = $data2[1] . $data2[2]; # combined TIR sequences of second member of pair
			my $side_sequences2 = $data2[0] . $data2[3]; # combined side sequences of second member of pair
			my $id2 = $data2[5]; # unique id of second member of pair
			my $tsd2 = $data2[8]; # TSD sequence of second member of pair
			my $te_sequence2 = ">" . $data2[4] . "#" . $id2 . "#" . $element . "\n" . $data2[6]; # fasta formated genomic sequence of second member of pair

			## if the two members of the pair are similar enough to be called the same

			my $TIR_IDENTITY_THRESHOLD = 0.1; # maximum p distance to call two tirs identical
			my $SIDE_IDENTITY_THRESHOLD = 0.7; # minimum distance between the side sequences to call them different 

			if ($location1 ne $location2) { # only continue if the two sequences are not just reverse complements of each other

				my $passed_tir_test = 0; # boolean, 0 if tirs are identical enough 1 if not
				my $passed_side_test = 0; # boolean, 0 if side sequence are the same 1 if they are different

				if (pdistance ($tirs1, $tirs2) < $TIR_IDENTITY_THRESHOLD) {
					$passed_tir_test = 1;
				}
				if (pdistance ($side_sequences1, $side_sequences2) > $SIDE_IDENTITY_THRESHOLD) {
					$passed_side_test = 1;
				}

				if ($passed_tir_test and $passed_side_test) {
					my $insertion_data1 = "$data[0]\t$data[7]\t$data[1]\t$data[2]\t$data[8]\t$data[3]\t$location1\t$te_sequence1";
					my $insertion_data2 = "$data2[0]\t$data2[7]\t$data2[1]\t$data2[2]\t$data2[8]\t$data2[3]\t$location2\t$te_sequence2";
					$insert_data{$insertion_data1} = "$data[1]\t$data[2]";
					$insert_data{$insertion_data2} = "$data2[1]\t$data2[2]";
				}
			}			
		}
	}

	## Output the data if any
	if (%insert_data) {
		print OUTPUT "\n";
		foreach my $key (sort {$insert_data{$a} cmp $insert_data{$b} } keys %insert_data) {
			my @data = split("\t", $key);
			my ($hitname, $evalue, $overlap) = blastseq($data[7]);
			my @data2 = split("\n", $data[7]);
			my $element_length = length($data2[1]);
			my $tir_length = tir_length($data2[1]);
			my $orf_lengths = size_longest_orfs($data2[1]);

			print OUTPUT "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$data[5]\t$hitname\t$evalue\t$overlap\t$element_length\t$tir_length\t$orf_lengths\t$data2[0]\t$data2[1]\n"; #all but location

		}
	}
}

# takes a sequence as input, does a blast search against a known db of TE proteins returns data about the hit 
sub blastseq {
	my ($sequence) = @_;

	# setup file for blast
	my $sequence_file = File::Temp->new( UNLINK => 1, SUFFIX => '.fa' ); # current sequence file
	my $blast_file = File::Temp->new( UNLINK => 1, SUFFIX => '.blast' ); # current sequence file
	open (SEQ, ">$sequence_file") or die "cannot open temporary file $sequence_file\n";
	print SEQ "$sequence\n";
	close SEQ;

	# execute the blast
	`fastx36 -E 1 $sequence_file $config{dbte} > $blast_file`;

	# interpret the blast results
	my $tename;
	my $te_length;
	my $evalue;
	my $overlap;

	if ($blast_file) {
		open (BLAST, $blast_file) or die "cannot open file $blast_file\n";
		while (my $line = <BLAST>) {
			if ($line =~ /^>>(\S+)\s+\((\d+)\saa\)/) {
				$tename = $1;
				$te_length = $2;
				$line = <BLAST>;
				my @data = split (" ", $line);
				$evalue = $data[-1];
				$line = <BLAST>;
				if ($line =~ /in\s+(\d+)\s+aa\s+overlap/) {
					$overlap = $1;
				}
				else {
					$overlap = "0";
				}
			}
			next if ($tename); # exit the loop once a hit has been found
		}
	}

	# return the results
	if ($tename) {
		my $percent_overlap = int(($overlap/$te_length) * 100);
		return ($tename, $evalue, $percent_overlap); 
	}
	else {
		return ("no_match", "N/A", "N/A");
	}
}

sub printUsage{

print STDOUT "DESCRIPTION: This program takes bed formated files from RepeatMasker and from genome TIR analysis and returns report\n";
print STDOUT "USAGE : combine-tir -r \"RepeatMasker .bed file\" -t \"TIR analysis .bed file\" -g \"Genome fasta file\" -d \"datbase of TE proteins\" -o \"output file\"
Options : 
    -r | rmbed		RepeatMasker output from RepeatExplorer file, bed formated (Mandatory)
    -t | tirbed		Whole genome TIR analysis, bed formated (Mandatory)
    -g | genome		Fasta formated genome (Mandatory)
    -d | dbte		Database of TE proteins (default $config{dbte})
    -o | out   		Name of ouput file (default \"out\")\n";
    exit;
}

#load a genome into a hash
sub uc_genometohash {
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
				$genome{$title} = uc($seq);
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
	$genome{$title} = uc($seq);

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

#take a te sequence file and returns length of TIRs
sub tir_length {
	(my $seq) = @_;
	
	my $MAX_TIR_LENGTH = 30; # maximum tir length
	my $MAX_MISSMATCH = 2; # maximum number of mismatches

 	# check that the sequence is long enough	
	if ((length $seq) < (2*$MAX_TIR_LENGTH)) {
		return ("too short");
	}
	
	my $miss = 0; # number of mismatches

	my $i = 0;
	while (($i <= $MAX_TIR_LENGTH) and ($miss < ($MAX_MISSMATCH+1))) {
		my $left_nuc = substr($seq, $i, 1);
		my $right_nuc = substr($seq, -$i - 1, 1);
		unless ($left_nuc eq rc($right_nuc)) {
			$miss++;
		}
		$i++;
	}

	$i--; # this is to account for that last $i++

	if ($i >= $MAX_TIR_LENGTH) {
		return(">=$i");
	}
	else {
		return($i);
	}
}

#get the longest ORFs
sub size_longest_orfs {
	(my $seq) = @_;

	my $NUMBER_OF_ORFs = 3; # number of longest orfs 

	# put the sequence into a temp file	
	my $orf_input_file = File::Temp->new( UNLINK => 1, SUFFIX => '.fa' );
	open (ORF, ">$orf_input_file") or die "cannot open file for writting $orf_input_file\n";
	print ORF ">temp\n$seq";

	# run the emboss getorf program
	my $orf_output_file = File::Temp->new( UNLINK => 1, SUFFIX => '.fa' );
	`getorf -sequence $orf_input_file -outseq $orf_output_file`;

	# interpret the data in the output file
	my @orflengths;
	my %orfs = uc_genometohash($orf_output_file);
	foreach my $orf_title (keys %orfs) {
		$orfs{$orf_title} =~ s/X//g;
		push @orflengths, (length $orfs{$orf_title});	
	} 

	#create the return text
	my @sorted_length = sort {$b <=> $a} @orflengths;
	my $return_text;
	for (my $i=0; $i < $NUMBER_OF_ORFs; $i++) {
		$return_text .= $sorted_length[$i] . " /";
	}
	chop $return_text;
	chop $return_text;

	return ($return_text);
	
}
exit;

