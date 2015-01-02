#!/usr/bin/perl

# Dec 2014 this script takes potential elements and searches them against a database of known elements

use strict;
use File::Temp();
use Getopt::Long;

##> Get input parameters 
my %config;
GetOptions(\%config,
	'tedb=s', # TE database
	'potele=s', # potential elements
	'genome=s', # genome sequence
	'dbgen=s', # genome blast db
	'out=s',
);

##> Check if no mandatory parameter is missing and set defaults
if (!exists $config{tedb})         {printUsage();}
if (!exists $config{potele})         {printUsage();}
if (!exists $config{genome})         {printUsage();}
if (!exists $config{dbgen})         {printUsage();}
if (!exists $config{out}) {$config{out} = "out";}
open (OUTPUT, ">$config{out}") or die "cannot open output file $config{out}\n";

# use fasta on the input file and store list of hits in %element_names
my $search_file = File::Temp->new( UNLINK => 1, SUFFIX => '.fst' ); # fasta output
my $search_result = `fasty36 -E 0.5 -T 12 -b 1 -m 8 $config{potele} $config{tedb} > $search_file`;
my %element_names;
open (INPUT, $search_file) or die;
while (my $line = <INPUT>) {
	my @data = split " ", $line;
	my $evalue = @data[-2];
	my $dbhitname = @data[1];
	$element_names{$data[0]} = "$dbhitname\t$evalue";
}
close INPUT;

# create preliminary report
my $filename = "prelim-report.xls";
open (PRELIMOUT, ">$filename") or die "cannot open file $filename for writting\n";
print PRELIMOUT "Element name\tHit name\tevalue\tsequence\n";
my %genomeseq = genometohash($config{genome});
my $blast_output = File::Temp->new( UNLINK => 1, SUFFIX => '.blst' ); # blast output
`mkdir motif_files`;
`blastn -outfmt 6 -num_alignments 1 -db $config{dbgen} -query $config{potele} -out $blast_output`;
open (INPUT, $blast_output) or die "cannot open file $blast_output\n";
my %positions;
while (my $line = <INPUT>) {
	my @data = split(" ", $line);
	my $name = $data[0];
	my $scaffold = $data[1];
	my $b1 = $data[8];
	my $b2 = $data[9];
	unless (exists $positions{$name}) {
		$positions{$name} = "$scaffold\t$b1\t$b2"; #holds the positions of the hits on the genome
	}
}
foreach my $element (keys %positions) {

	# getting the sequence from the genome
	my $sequence;
	my @data = split(" ", $positions{$element});
	if ($data[1] < $data[2]) {
		$sequence = uc(substr($genomeseq{$data[0]}, $data[1], ($data[2] - $data[1])));
	}
	if ($data[2] < $data[1]) {
		$sequence = substr($genomeseq{$data[0]}, $data[2], ($data[1] - $data[2]));
		$sequence = uc(rc($sequence));
	}

	# finding TE motifs
	my @data3 = split("#", $element);
	my $number_of_motifs = findmotifs($data3[0], $sequence);

	my @data2 = split (" ", $element_names{$element});
	print PRELIMOUT "$element\t$data2[0]\t$data2[1]\t$sequence\n";
}
close PRELIMOUT;

#print the those sequences that matched
my %proposed_elements = genometohash($config{potele});
foreach my $key (keys %element_names) {
print "ha\n";
	print OUTPUT ">$key\n";
	print OUTPUT "$proposed_elements{$key}\n";
}

close OUTPUT;

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

sub printUsage{

print STDOUT "DESCRIPTION: This program takes a FASTA formated file of potential element and compares it to a database of known elements\n";
print STDOUT "USAGE : perl search_potential_elements.pl -t \"known TE database, fasta file\" -p \"potential elements, fasta file\" -g \"genome sequence\" -o \"output file\"
Options : 
    -t | tedb		Database of known elements, fasta formated (Mandatory)
    -p | potele		List of potential elements, fasta formated (Mandatory)
    -g | genome		genome sequence, fasta formated (Mandatory)
    -d | dbgen		genome blast database (Mandatory)
    -o | out   		Name of ouput file (default \"out\")\n";
    exit;
}

sub rc {
    my ($sequence) = @_;
    $sequence = reverse $sequence;
    $sequence =~ tr/ACGTRYMKSWacgtrymksw/TGCAYRKMWStgcayrkmws/;
    return ($sequence);
}

sub findmotifs {
	my ($name, $sequence) = @_;
	my $pfamREPET="/home/arensburger/db/repet/ProfilesBankForREPET_Pfam26.0_GypsyDB.hmm";
	
	my $sequence_file = File::Temp->new( UNLINK => 1, SUFFIX => '.fas' ); # contains sequence to translate
	my $translated_file = File::Temp->new( UNLINK => 1, SUFFIX => '.fas' ); # translated file
	open (OUTPUT2, ">$sequence_file") or die;
	print OUTPUT2 ">temp\n";
	print OUTPUT2 "$sequence\n";
	`transeq $sequence_file $translated_file -frame=6 2>&1`; # translate the DNA into 6 frames
	my $hmm_output = File::Temp->new( UNLINK => 1, SUFFIX => '.fas' ); # translated file
	`hmmscan --tblout $hmm_output $pfamREPET $translated_file`;
	my $number_of_hits =  `grep -v "^#" $hmm_output -c`;
	chomp $number_of_hits;
	if ($number_of_hits) {
		my $output_file = $name . "-hmmprofile.txt\n";
		`cp $hmm_output motif_files/$output_file`;
	}
	return($number_of_hits);
}

