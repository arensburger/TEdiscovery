#!/usr/bin/perl
# Given a genomic location find possible TSD and TIR combination

use strict;
use Bio::SearchIO;
use Bio::Index::Fasta;
use lib "/home/peter/scripts/elesearch";
use agutils;
use Getopt::Long;

#constants to be modified for every new data set
my $MAX_ELEMENT_SIZE = 10000; #supercedes the $SEARCH_DISTANCE parameter
my $SEARCH_DISTANCE = 5000; #nucleotides away from hit to search for TSDs and TIRs
my $TIR_LEN = 11;
my $TIR_MISSMATCH = 2;
my $TSD_MISSMATCH = 2;
my $TSD_LEN = 2; #length of the TSDs
my $searchseq; #nucleotide sequence with a blast hit to an element, search this for TSDs and TIRs 
my $eleseq; #sequence of the potential element
my $seq_name;
my $c1;
my $c2;
my $GENOME = "/data/genomes/Analbi/anop_albi_albi9_a.i1.scaffolds.fasta";
my $INPUT_POSITIONS; #file that holds the positions of the input sequences
my $note; #name to give output file

# setup to modify constants on the command line
GetOptions(
	'i:s'	=> \$INPUT_POSITIONS,
	't:s'	=> \$TSD_LEN,
	'g:s'	=> \$GENOME,
);
#check the inputs
unless ($INPUT_POSITIONS) {
	die "USAGE: perl find_tirs.pl -i <REQUIRED: file with input positions> -t <OPTIONAL: length of TSD> -g <OPTIONAL: genome file>\n";
}

#read input file
open (INPUT, $INPUT_POSITIONS) or die "cannot open file $INPUT_POSITIONS\n";
while (my $line = <INPUT>) {
	my %output_data; #holds the sequence length as key and output line as value, used for sorting the output
	if ($line =~ /^(\S+):(\d+)-(\d+)\t(\S+)\s/) {
		$seq_name = $1;
		$c1 = $2;
		$c2 = $3;
		$note = $4;

		if ($c1 > $c2) {
			my $temp = $c1;
			$c1 = $c2;
			$c2 = $temp;
		}
	}
	else {
		die "cannot read input line\n$line";
	}

#	open (OUTPUT, ">>$seq_name:$c1-$c2-tsdtir.txt") or die "cannot make output file\n";
       open (OUTPUT, ">>$note") or die "cannot make output file\n";

	#calculate search distance based on maximum element length then get boundaries
	$SEARCH_DISTANCE = int(($MAX_ELEMENT_SIZE - ($c2 - $c1))/2);
	my $b1 = $c1 - $SEARCH_DISTANCE;
	my $b2 = $c2 + $SEARCH_DISTANCE;
	if ($b1 < 1) {
		$b1 = 1;
	}
	(my $title, my $searchseq) = ext_genome($GENOME, $seq_name, $b1, $b2);
	$searchseq = uc($searchseq);

	my $left_seq = substr ($searchseq, 1, ($c1 - $b1)); #sequence to look for TSDs and TIRs on left sid
	my $right_seq = substr ($searchseq, (length $left_seq) + ($c2 - $c1), $SEARCH_DISTANCE);  #same on the right

	#look for all possible TIR matches
	for (my $i=0; $i<(length $left_seq) - $TIR_LEN; $i++) {
		my $ltir = substr($left_seq, $i-1, $TIR_LEN);
		my $ltsd = substr($left_seq, $i-$TSD_LEN - 1, $TSD_LEN);
		my $lpos = $i - $TSD_LEN;
		unless ($ltsd =~ /n|N/) {		
			for (my $j=0; $j<(length $right_seq) - $TIR_LEN; $j++) {
				my $rtir = substr($right_seq, $j, $TIR_LEN);
				my $rtsd = substr($right_seq, $j+$TIR_LEN, $TSD_LEN);
				my $rpos = (length $searchseq) - (length $right_seq) + $j + $TIR_LEN + $TSD_LEN;
				if (itrmatch_n($ltsd, rc($rtsd), $TSD_MISSMATCH) == 1) { #if the TSD match
					my $continue_processing = 1; #boolean, used to stop looking at this element
					if ($TSD_LEN == 2) { #if we're looking at 2bp TSD only continue if we have TA TSDs
						unless (($ltsd eq "TA") and ($rtsd eq "TA")) {
							$continue_processing = 0;
						}
					}
					if ((itrmatch_n($ltir, $rtir, $TIR_MISSMATCH) == 1) and ($continue_processing)) { #the TIRs match
						my $rc_rtir = rc $rtir;
						my $seq = substr($searchseq, $lpos, ($rpos - $lpos - 1));
						my $corrected_lpos = $b1 + $lpos;
						my $corrected_rpos = $corrected_lpos + length ($seq);	
		
						#calculate TSD match
						my $tsd_match;
						for (my $l=0; $l < length $ltsd; $l++) {
							my $t1 = substr($ltsd, $l, 1);
							my $t2 = substr($rtsd, $l, 1);
							if ($t1 eq $t2) {
								$tsd_match++;
							}
						}
		
						#try to see how long the TIR still work
						my $exit_loop = 0; #boolean exit once the TIRs don't match anymore
						my $k=0;
						my $max_tir_len = $TIR_LEN;
						my $long_ltir;
						my $long_rtir;
						while ($exit_loop == 0) {
							$k++;
							$long_ltir = substr($left_seq, $i-1, ($TIR_LEN + $k));
							$long_rtir = substr($right_seq, ($j-$k), ($TIR_LEN + $k));
							if (itrmatch_n($long_ltir, $long_rtir, $TIR_MISSMATCH) == 1) {
								$max_tir_len++;
							}
							else {
								$exit_loop = 1;
							}
							if ($k > 30) {
								$exit_loop = 1;
							}
						}
	
						#print results
						my $len_seq = length $seq;
						$output_data{$len_seq} = "$seq_name:$corrected_lpos-$corrected_rpos\t$len_seq\t$tsd_match\t$max_tir_len\t$ltsd\t$rtsd\t$long_ltir\t$long_rtir\n";
					}
				}
			}
		}
	}

	#print output sorting by size
	for my $key ( sort {$a<=>$b} keys %output_data) {
           	print OUTPUT "$output_data{$key}";
	}
	close OUTPUT;
}
close INPUT;
