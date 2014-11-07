#!/usr/bin/perl
use strict;
use Getopt::Long;
use agutils;
use File::Temp();

my $EXPAND_DISTANCE = 1000; # bp on each side of the rm hit to expand
my $MIN_AA_OVERLAP = 20; # minimam number of aa for overlap to be real
my $MAX_PRINT_COUNT = 1;
my $TE_blast_db = "./invrepaa/invrep-aa.fa"; # name of blast database (as used for -db blast argument) with TEs
my $blastoutname = "blast-files3";

my $rm_filename; #reapeatmasker .out file name
my $genome_filename;

# read and check the inputs
GetOptions(
	'i:s'	=> \$rm_filename,
	'g:s'	=> \$genome_filename,
);
unless (($rm_filename) and ($genome_filename)) {
	die "USAGE: perl rm-te-parse -i <REQUIRED: RepeatMasker .out file> -g <REQUIRED: genome fasta file>\n";
}

# load the reference genome into memory
print STDERR "loading genome into memory...\n";
my %genome = genometohash($genome_filename);

# read the RepeatMasker file
my %seqs; #tename as key as value array with [0] location on 
print STDERR "reading rm file...\n";
open (INPUT, $rm_filename) or die "$rm_filename: $!\n";
`rm -r $blastoutname`;
`mkdir $blastoutname`;

my %allseq; # holds all the summary data, with sequence name as key
my %nucseq; # holds the sequences associated with each RM sequencs as key

# figure out RepeatMasker file size
(my $rmsize) = `wc -l $rm_filename` =~ /\A([^:\s]+)/;
my $i=0; # line counter
<INPUT>;
<INPUT>;
<INPUT>;
while (my $line = <INPUT>) {
	# read the rm data
	my @data = split(" ", $line);
	my $contig = $data[4];
	my $b1 = $data[5];
	my $b2 = $data[6];
	my $left = $data[7];
	$left = substr($left, 1, -1); 
	my $ori = $data[8];
	my $seqname = $data[9];
	my $class = $data[10];

	# figure out the bounds to extract from the contig
	my $contigb1 = $b1 - $EXPAND_DISTANCE;
	if ($contigb1 < 1) { 
		$contigb1 = 1;
	}
	my $contigb2;
	if ($left < $EXPAND_DISTANCE) { 
		$contigb2 = $b2 + $left;
	}
	else {
		$contigb2 = $b2 + $EXPAND_DISTANCE;
	}
	# extract the sequence if it's a class II sequence
	if ($class =~ /^DNA/) {
               my $sequence = substr($genome{$contig}, $contigb1, ($contigb2 - $contigb1));
		# set up and run the blasts
		my $blast_file = File::Temp->new( UNLINK => 1, SUFFIX => '.tbl' ); # blast output
		my $sequence_file = File::Temp->new( UNLINK => 1, SUFFIX => '.tbl' ); # current sequence file
		my $outputdir = "$blastoutname/$seqname";
		
		open (OUTPUT, ">$sequence_file") or die "cannot open output file $sequence_file\n";
		print OUTPUT ">$contig:$contigb1-$contigb2\n$sequence\n";
		close OUTPUT;
		
		`fastx36 -E 1 $sequence_file $TE_blast_db > $blast_file`;
		open (INPUT2, $blast_file) or die "cannot open file $blast_file\n";
		my $print_count = 0; # keeps track of how many blast hits have been outputed 
		while ((my $line2 = <INPUT2>) and ($print_count < $MAX_PRINT_COUNT)) {
			my $tename;
			my $te_length;
			my $evalue;
			my $overlap;
			if ($line2 =~ /^>>(\S+)\s+\((\d+)\saa\)/) {
				$tename = $1;
				$te_length = $2;
				$line2 = <INPUT2>;
				my @data = split (" ", $line2);
				$evalue = $data[-1];
				$line2 = <INPUT2>;
				if ($line2 =~ /in\s+(\d+)\s+aa\s+overlap/) {
					$overlap = $1;
				}
				else {
					$overlap = "0";
				}
				$print_count = 1;

                        	if ($overlap >= $MIN_AA_OVERLAP) {
					unless (-e $outputdir) {
						`mkdir $outputdir`;
					}
					my $name = "fastx36" . "-" . $contig . ".txt";
			                `cp $blast_file $outputdir/$name`;
					`cat $sequence_file >> $blastoutname/genome_sequences.fa`; 
					my $percent = $overlap/$te_length;

					$allseq{$seqname} .= "$seqname\t$contig:$contigb1-$contigb2\t$ori\t$tename\t$evalue\t$overlap\t$te_length\t$percent\n";	
					$nucseq{$seqname} .= ">$contig:$contigb1-$contigb2\n$sequence\n";	
                                	# print "$seqname\t$contig\t$contigb1\t$contigb2\t$ori\t$tename\t$evalue\t$overlap\t$te_length\t$percent\n";
                        	}
			}

		} 
		close INPUT2;
	}
}

print "RM input name\tgenome location\tdb hit name\tevalue\toverlap\tdb TE length\tratio\n";
foreach my $seq (keys %allseq) {
	my $outputdirname = $blastoutname . "/" . $seq;
	print "$allseq{$seq}\n";
	my $dotterseq = "$outputdirname/dotter-seq.fa";
	my $dotteroutput = "$outputdirname/dotter-output";
	open (DOTOUT, ">$dotterseq") or die "cannot create file $dotterseq\n";
	print DOTOUT "$nucseq{$seq}";
	close DOTOUT;
	`dotter_ul1 -b $dotteroutput $dotterseq $dotterseq`;
}
exit;
