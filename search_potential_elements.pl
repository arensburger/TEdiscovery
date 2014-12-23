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
	'out=s',
);

##> Check if no mandatory parameter is missing and set defaults
if (!exists $config{tedb})         {printUsage();}
if (!exists $config{potele})         {printUsage();}
if (!exists $config{out}) {$config{out} = "out";}
open (OUTPUT, ">$config{out}") or die "cannot open output file $config{out}\n";


#my $TE_database = $ARGV[0]; # fasta file with known TE amino acid sequences
#my $possible_element_file = $ARGV[1]; #fasta file of proposed elements

# create temporary blast database
#my $tedb_dir = File::Temp->newdir(); # TE database directory
#`makeblastdb -in $config{tedb} -dbtype prot -out $tedb_dir/TEdb`;

# use fasta on the input file and store list of hits in %element_names
my $search_file = File::Temp->new( UNLINK => 1, SUFFIX => '.fst' ); # fasta output
my $search_result = `fasty36 -E 0.5 -T 12 -b 1 -m 8 $config{potele} $config{tedb} > $search_file`;
my $output_name = $config{potele} . ".fasty";
`cp $search_file $output_name`;
my %element_names;
open (INPUT, $search_file) or die;
while (my $line = <INPUT>) {
	my @data = split " ", $line;
	$element_names{$data[0]} = 0;
}

#print the those sequences that matched
my %proposed_elements = genometohash($config{potele});
foreach my $key (keys %element_names) {
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
print STDOUT "USAGE : perl search_potential_elements.pl -t \"known TE database, fasta file\" -p \"potential elements, fasta file\" -o \"output file\"
Options : 
    -t | tedb		Database of known elements, fasta formated (Mandatory)
    -p | potele		List of potential elements, fasta formated (Mandatory)
    -o | out   		Name of ouput file (default \"out\")\n";
    exit;
}
