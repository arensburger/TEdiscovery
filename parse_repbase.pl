#!/usr/bin/perl

# Nov 2014.  Takes Repbase EMBL file and returns part

open (INPUT, $ARGV[0]) or die "cannot open file $ARGV[0]\n";

while (my $line = <INPUT>) {
	if ($line =~ /^FT\s+\/product="(\S+)"/) {
		my $name = $1;
		my $seq;
		do {	
			$line = <INPUT>;
			if ($line =~ /\/translation="(\S+)/) {
				$seq .= $1;
			}
			elsif ($line =~ /FT\s+\//) {
			}
			else {
				if ($line =~ /^FT\s+([A-Z]+)/) {
					$seq .= $1;
				}
			} 
		} until ($line =~ /^XX/);
		print ">$name\n$seq\n";
	}
	
}
