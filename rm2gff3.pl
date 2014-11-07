#!/usr/bin/perl
# April 2014 Takes output of RE blast the output

use strict;
use File::Temp ();
use Getopt::Long;

my %config;
GetOptions(\%config,
	'rmout=s',
	'out=s',
);

##> Check if no mandatory parameter is missing and set defaults
if (!exists $config{rmout})         {printUsage();}
if (!exists $config{out}) {$config{out} = "out";}

#### Main program ###
open (RM, "$config{rmout}") or die "cannot open file $config{rmout}\n";
open (OUT, ">$config{out}") or die "cannot open output file $config{out}\n";
<RM>;
<RM>;
<RM>;
#print OUT "te_name\tgene_type\tgene_creation_date\tgene_comments\tte_owner\tte_scaffold\tte_start\tte_end\tte_strand\tte_superfamily\tte_family\tfull_or_partial\n";
while (my $line = <RM>) {
	my @data = split (" ", $line);
	my $scaffold = $data[4]; 
	my $b1 = $data[5]; 
	my $b2 = $data[6];
	my $ori = $data[8];
	my $name = $data[9];
	if ($ori eq "C") {
		$ori = "-";
	}
	my $ele = $data[10];
	my $sf;
	my $fam;	
	if ($ele =~ /(\S+)\/(\S+)/) {
		$sf = $1;
		$fam = $2
	}
	else {
		$sf = $ele;
		$fam = "unknown"
	}

	unless ( ($sf eq "Low_complexity") or ($sf eq "Simple_repeat") or ($sf eq "Satellite")) {
#	if ($sf) {
#		print OUT "$name\tTransposon\t2014-10-16\tpreliminary annotation\tparensburger\t$scaffold\t$b1\t$b2\t$ori\t$sf\t$fam\n";
		my $name2 = "$name" . "#" . "$sf" . "#" . "$fam";
		print OUT "$scaffold\t$b1\t$b2\t$name2\t0\t$ori\n";
	}
}
close RM;
close OUT;


sub printUsage{

print STDOUT "DESCRIPTION: This program takes the .out file of a RepeatMakser run and converts it into a gff3 file\n";
print STDOUT "USAGE : rm2gff3.pl -r \"RepeatMasker .out file\" -o \"output file\"
Options : 
    -r | rmout		RepeatMasker .out file (Mandatory)
    -o | out   		Name of ouput file (default \"out\")\n";
    exit;
}
exit;
