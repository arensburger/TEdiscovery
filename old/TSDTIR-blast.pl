#!/usr/bin/perl
# April 2014 Takes output of RE blast the output

use strict;
use File::Temp ();
use Bio::SearchIO;
use Getopt::Long;

my $NUMBEROFCPUS = 8;
my $REPBASEDATABASE = "/home/parensburge/db/Repbase_tefam/reptef"; # formated for use with blastn

### read and check the inputs
my $infile; # input fasta file
my $outputdirectory; 
my $bedfile; # bed file that holds all the tir ends
GetOptions(
	'in:s'   => \$infile,
	'b:s'	=> \$bedfile,
	'o:s'	=> \$outputdirectory,
);
unless ($infile and $outputdirectory and $bedfile) {
	die ("usage: perl TSDTIR-blast -in <input fasta file, REQUIRED> -b <bed file with TIR ends, REQUIRED> -o <output directory, REQUIRED>");
}
`mkdir -p $outputdirectory`;
if ($?){
	die;
} 

# get all the tir ends for each element
my %tir; # holds the sequence bounds as key and tirs (tab separated) as value
my %tirfreq; # holds the tir sequence as key and frequency (both orientations as value)
my %tirid; # holds tirs sequence as key and unique id as value
open (INPUT, $bedfile) or die "cannot open bed file $bedfile\n";
my $i=0;
while (my $line = <INPUT>) {
	if ($line =~ /^(\w+)\s(\w+)\s\w+\s\w+\s\w+\s(\w+)\s.+tir1=(\w+);tir2=(\w+)/) {
		my $loc = $1 . ":" . $2 . "-" . $3;
		my $tirseq = "$4\t$5";
		$tir{$loc} = $tirseq;
		# update te frequency
		if (exists $tirfreq{rc ($tirseq)}) {
			$tirfreq{rc ($tirseq)} += 1;
		}	
		else {
			$tirfreq{$tirseq} += 1;
		}

		# update te id
		unless ((exists $tirid{$tirseq}) or (exists $tirid{rc ($tirseq)})){
			$tirid{$tirseq} = $i;
			$i++;
		}
		
	}
	else {
		warn "cannot read line in bed file\n$line";
	}
}
close INPUT;

## perform the blast using nucleotide blast
#my $rm_blastoutput = "$outputdirectory/blastn_output.txt";
#`blastn -query $infile -db $REPBASEDATABASE -out $rm_blastoutput -outfmt 1 -evalue 0.01 -num_alignments 1 -num_threads $NUMBEROFCPUS 2>/dev/null`;

# get all the orfs and blast them
my %tirorf; # tir location as key and list of elements that match as value
my $getorfinfile = File::Temp->new( UNLINK => 1, SUFFIX => '.fa' );
my $orfout = File::Temp->new( UNLINK => 1, SUFFIX => '.fa' );
my $orfblastoutput = "$outputdirectory/orf_tblastn.txt";
# get orf does not deal with colons, need to replace them with underscore
open (INPUT, $infile) or die "cannot open file\n";
open (OUTPUT, ">$orfblastoutput") or die;
while (my $line = <INPUT>) {
	$line =~ s/:/_/g;
	print OUTPUT $line;
}
close INPUT;
close OUTPUT;
`getorf -minsize 300 -sequence $orfblastoutput -outseq $orfout`;
`tblastn -query $orfout -db $REPBASEDATABASE -out $orfblastoutput -num_threads $NUMBEROFCPUS -outfmt 6 -num_alignments 1 -evalue 0.001`;
open (INPUT, $orfblastoutput) or die "cannot open file $orfblastoutput\n";
while (my $line = <INPUT>) {
	my @data = split(" ", $line);

	# recreate the original title
	my $title = $data[0];
	my $new_title;
	if ($title =~ /(\S+)_(\d+-\d+)_/){
		$new_title = $1 . ":" . $2;
	}
	else {
		die "cannot read line $title";
	}

	$tirorf{$new_title} = $data[1];
}

# perform the blast on RM using tblastx
my $tblastxoutput = "$outputdirectory/tblastx.xml";
`tblastx -query $infile -db $REPBASEDATABASE -out $tblastxoutput -outfmt 5 -evalue 0.1 -num_alignments 1 -num_threads $NUMBEROFCPUS `;
#2>/dev/null
my $tblastxsimple = "$outputdirectory/tblastx.xls";
open (OUTPUT, ">$tblastxsimple") or die "cannot create file $tblastxsimple\n";
my $searchin = new Bio::SearchIO( -tempfile => 1,
		 		 -format => 'blastxml',
		 		 -file   => $tblastxoutput);
while( my $result = $searchin->next_result ) {
	my $query_length = $result->query_length;
	my $query_name = $result->query_description;
	while (my $hit = $result->next_hit) {
		my $hit_name = $hit->name;
		my $hsp = $hit->next_hsp;
		my $evalue = $hsp->evalue;
		my $description = $hit->description;
		my $start_query = $hsp->start('query');
		my $end_query = $hsp->end('query');
		my $proportion = ($end_query - $start_query)/$query_length;
		print OUTPUT "$query_name\t$hit_name\t$evalue\t$query_length\t$start_query\t$end_query\t$proportion\t$tir{$query_name}\t$tirfreq{$tir{$query_name}}\t$tirid{$tir{$query_name}}\t$tirorf{$query_name}\n";

	}
}
close OUTPUT;

#### SUBROUTINES ########
#reverse complement
sub rc {
    my ($sequence) = @_;
    $sequence = reverse $sequence;
    $sequence =~ tr/ACGTRYMKSWacgtrymksw/TGCAYRKMWStgcayrkmws/;
    return ($sequence);
}
