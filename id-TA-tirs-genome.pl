#!/usr/bin/perl
# Takes a fasta file and id's TIRs

use strict;

my $TSD = "TA";
my $TIR_LENGTH = 11;
my $MAX_ELE_LENGTH = 2500;
my $WINDOW_LENGTH = 5000;
my $MIN_ELE_LENGTH = 500;
my $MAX_TIR_MISSMATCH = 3;
my $MIN_PERCENT_GC_TIR = 25; # mininum percent GC in TIRs

my %hitloc; #holds the location of potential TEs as key and differences between the TEs as value

# load the sequence
my %sequence = genometohash($ARGV[0]); #converts to upper case too

foreach my $seqname (keys %sequence) {
	# identify all the TSDs in this sequence and populate the arrays
	my @fpos; #holds the position of the first bp all forward tirs
	my @fseq; #holds the sequence of all the forward tirs
	my @rpos; #same as above
	my @rseq;
	print STDERR "$seqname\n";	
	#identify all the TSDs and associated TIR sequences
	while ($sequence{$seqname} =~ /$TSD/g) {
		my $position = pos($sequence{$seqname}) - 1;
		push(@fpos, $position);
		push (@fseq, substr($sequence{$seqname}, $position + 1, $TIR_LENGTH));
		push (@rseq, substr($sequence{$seqname}, $position - $TIR_LENGTH - 1, $TIR_LENGTH));
	}

	#identify and compare all thing in a window
	for (my $i=0; $i < length ($sequence{$seqname}); $i=$i+$MAX_ELE_LENGTH) {
		my @local_fpos=(); #holds the array index of all fpos in current window
		my @local_fseq=(); #holds the array index of all fpos in current window
		my @local_rseq=(); #holds the array index of all fpos in current window
		#identify all local positions for this window
		my $notpassed = 1; #boolean true until index is higher than window
		my $j=0;
		while ($notpassed) {
			if (($fpos[$j] >= $i) and ($fpos[$j] <= ($i + $WINDOW_LENGTH))) {
				push (@local_fpos, $fpos[$j]);
				push (@local_fseq, $fseq[$j]);
				push (@local_rseq, $rseq[$j]);
			}
			if ($fpos[$j] > ($i + $WINDOW_LENGTH)) {
				$notpassed = 0;
			}
			if ($j >= scalar @fpos) {
				$notpassed = 0;
			}
			$j++;
		}

		#compare all the forward to reverse sequences
		for (my $k=0; $k<scalar @local_fpos; $k++) { #scrolling through forward sequences
	
			# test for too AT rich
			my $f_gc = ($local_fseq[$k] =~ tr/C//); 
			$f_gc += ($local_fseq[$k] =~ tr/G//);
			my $f_gc_percent = ($f_gc/$TIR_LENGTH) * 100;
			next if ($f_gc_percent < $MIN_PERCENT_GC_TIR);

			for (my $l=$k+1; $l < scalar @local_fpos; $l++) {
				my $seqsize = $local_fpos[$l] - $local_fpos[$k];
				next if ($seqsize < $MIN_ELE_LENGTH); #exit if element is too small

				# test for too AT rich
				my $r_gc = ($local_rseq[$l] =~ tr/C//); 
				$r_gc += ($local_rseq[$l] =~ tr/G//);
				my $r_gc_percent = ($r_gc/$TIR_LENGTH) * 100;
				next if ($r_gc_percent < $MIN_PERCENT_GC_TIR);

				#count the differences
				my $differences = 0; #number of differences
				my $rc_local_rseq = rc($local_rseq[$l]); #rc of reverse sequence
				for (my $m=0; $m<$TIR_LENGTH; $m++) {
					if ((substr($local_fseq[$k], $m, 1)) ne (substr($rc_local_rseq, $m, 1))) {
						$differences++;
					}
				} 

				if ($differences < $MAX_TIR_MISSMATCH) {
					my $adjusted_left_pos = $local_fpos[$k] - 1;
					my $adjusted_right_pos = $local_fpos[$l] + length($TSD) - 1; # placing the right postion at the end of the TSD
					my $loc = "$seqname\t" . $adjusted_left_pos . "\t" . $adjusted_right_pos . "\t" . "TA_11TIR";
					$hitloc{$loc} = $differences;
#					print "$seqname\t$adjusted_left_pos\t$local_fseq[$k]\t$adjusted_right_pos\t$local_rseq[$l]\t$differences\n";	
				}
			}
			
		} 
	}
} 

#print results
foreach my $key (keys %hitloc) {
	print "$key\n";
}	




sub rc {
    my ($sequence) = @_;
    $sequence = reverse $sequence;
    $sequence =~ tr/ACGTRYMKSWacgtrymksw/TGCAYRKMWStgcayrkmws/;
    return ($sequence);
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
			$seq .= uc($line);
		}
	} 
	$genome{$title} = $seq;

	return (%genome);
}
