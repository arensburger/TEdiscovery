#!/usr/bin/perl

# Nov 2014 this script takes potential elements and blasts them against a database of known elements

use strict;
use File::Temp();

my $TE_database = $ARGV[0]; # fasta file with known TE amino acid sequences
my $possible_element_file = $ARGV[1]; #fasta file of proposed elements

# create temporary blast database
my $tedb_dir = File::Temp->newdir(); # TE database directory
`makeblastdb -in $TE_database -dbtype prot -out $tedb_dir/TEdb`;

# blast the input file
my $blast_file = File::Temp->new( UNLINK => 1, SUFFIX => '.tbl' ); # blast output
`blastx -db $tedb_dir/TEdb -query $possible_element_file -evalue 0.5 -out $blast_file -outfmt 6`;

#identify the elements that passed the test
my %element_names;
open (INPUT, $blast_file) or die;
while (my $line = <INPUT>) {
	my @data = split " ", $line;
	$element_names{$data[0]} = 0;
}

#print the those sequences that matched
my %proposed_elements = genometohash($possible_element_file);
foreach my $key (keys %element_names) {
	print ">$key\n";
	print "$proposed_elements{$key}\n";
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
