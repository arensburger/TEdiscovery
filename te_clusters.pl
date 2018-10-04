# oct 2010 this script replicates what Brennecke did in 2007 looking for small RNA clusters
# Feb 2013 takes a .sam file a input (presumable small RNAs mapped to a genome), returns a file with the following columns
# Col 0: unique cluster ID
# Col 1: contig id from sam file
# Col 2: bound 1 of cluster
# Col 3: bound 2 of cluter
# Col 4: cluster size (difference between Col 2 and 3)
# Col 5: number of piRNAs on the positive strand
# Col 6: number of piRNAs on the negative strand
# Col 7: number of unique piRNAs that could be generated from this cluster
 
use strict;
use Getopt::Long;

# declare some general variables
my %clusterplus; #holds the name of potential clusters as key and mapping frequency on plus strand as value
my %clusterminus; #holds the name of potential clusters as key and mapping frequency on plus strand as value
my %contnamesize; #holds the name of the contigs hit as key and the highest cluster as value, used to know how far too look on a supercontig
my %bounds; #holds as key the cluster name and as value [0] location of the lowest "location" RNA, [1] highest locations, used to get exact bound on the cluster location
my %seqcluster; #holds the sequence of the matches as key and as value the cluster it matches to 

#read and check the inputs values
my $WINDOWSIZE = 5000; #window size to look for small RNAs
my $inputname; #sam formated input file
my $outputname;
my $minlen = 24; #minimum length of piRNA
my $maxlen = 32; #maximum length of a piRNA
GetOptions(
	'w:s'   => \$WINDOWSIZE,
	'o:s'	=> \$outputname,
	'i:s'   => \$inputname,
	'l:s'	=> \$minlen, 
	'm:s'	=> \$maxlen
);

unless ($inputname) {
	die "usage: perl te_cluster.pl -i <REQUIRED: input file in .sam format>\n\t-o <OPTIONAL: output file; default stdout>\n\t-w <OPTIONAL: window size; default $WINDOWSIZE>\n\t-l <OPTIONAL: mininum length of piRNA; default $minlen>\n\t-m <OPTIONAL: maximum length of piRNA; default $maxlen>\n";
}
open (INPUT, $inputname) or die "cannot open file $inputname\n";
if ($outputname) {
	open (OUTPUT, ">$outputname") or die "cannot open output name $outputname\n";
}

#### start the analysis proper

#parse the sam input file
while (my $line = <INPUT>) {
   if ($line =~ /^\S+\s(\d+)\s(\S+)\s(\d+)\s\d+\s\S+\s\S+\s\S+\s\S+\s(\S+)\s/) {
    my $orientation = $1;
    my $contig = $2;
    my $location = $3;
    my $seq = $4;
    my $length = length $seq;
    if (($orientation != 4) && ($length >= $minlen) && ($length <= $maxlen)) { #continue only maps and size is ok
      #update the cluster density
      my $index = int($location/$WINDOWSIZE);
      my $clustername = $contig . "_" . $index;
      if ($orientation eq "0") {
      	$clusterplus{$clustername} += 1;
      }
      else {
	$clusterminus{$clustername} += 1;
      }
      unless ($index < $contnamesize{$contig}) {
	     $contnamesize{$contig} = $index;
      } 

      #update the highest and lowest location for this cluster
      if (exists $bounds{$clustername}) {
	if ($location < $bounds{$clustername}[0]) {
	  $bounds{$clustername}[0] = $location;
	}
	if ($location > $bounds{$clustername}[1]) {
	  $bounds{$clustername}[1] = $location;
	}
      }
      else {
	$bounds{$clustername}[0] = $location;
	$bounds{$clustername}[1] = $location;
      }

      $seqcluster{$seq} = $clustername;
    }
  }
  else {
#	warn "could not read line\n$line";
  }
}

#determine the number of different piRNAs mapped by each cluster
my %diffpiRNAs; #holds the cluster name as key and number of diffrent piRNA generated from that cluster as value
foreach my $key (keys %seqcluster) {
  $diffpiRNAs{$seqcluster{$key}} += 1;
}

#find cluster, go through each contig and group the indexes that have hits
my %cluster; #holds the number of the cluster as key and an array as value [0] contig of cluster, [1] start of cluster, [2] end of cluster, [3] total number of hits in the on + [4] total number of hits on - [5] number of piRNAs mapped by this cluster
my $clusnum=0; #current cluster number
foreach my $contig (keys %contnamesize) {
	my $inclus = 0; #boolean for in cluster or not
	for (my $i = 1; $i <= $contnamesize{$contig}; $i++) {
		my $location = $contig . "_" . $i;
		my $count = $clusterplus{$location} + $clusterminus{$location};
		if ($count > 0) { #this cluster is populated
			if ($inclus == 0) { #currently not adding to a cluster
				$clusnum++;
				$cluster{$clusnum}[0] = $contig;
				$cluster{$clusnum}[1] = $bounds{$location}[0];
				$cluster{$clusnum}[2] = $bounds{$location}[1];
				$cluster{$clusnum}[3] = $clusterplus{$location};
				$cluster{$clusnum}[4] = $clusterminus{$location};
				$cluster{$clusnum}[5] = $diffpiRNAs{$location};

				$inclus=1; #now we are in a cluster
			}
			else {
				$cluster{$clusnum}[2] = $bounds{$location}[1];
				$cluster{$clusnum}[3] += $clusterplus{$location};
				$cluster{$clusnum}[4] += $clusterminus{$location};
				$cluster{$clusnum}[5] += $diffpiRNAs{$location};

			}
		}
		elsif ($inclus == 1) { #we are adding to a cluster but now count is now 0
			$inclus = 0;
		}
	}
}

foreach my $key (keys %cluster) {
	my $clustersize = $cluster{$key}[2] - $cluster{$key}[1];

	#because I don't intialize $cluster{$clusnum}[3] and $cluster{$clusnum}[4] set these to zero if they are blank
	if (length($cluster{$key}[3]) == 0) {
		$cluster{$key}[3] = 0;
	}
	if (length($cluster{$key}[4]) == 0) {
		$cluster{$key}[4] = 0;
	}
#print "$cluster{$clusnum}[2]\n";

	if ($outputname) {
		print OUTPUT "$key\t$cluster{$key}[0]\t$cluster{$key}[1]\t$cluster{$key}[2]\t$clustersize\t$cluster{$key}[3]\t$cluster{$key}[4]\t$cluster{$key}[5]\n";
	}
	else {
		print "$key\t$cluster{$key}[0]\t$cluster{$key}[1]\t$cluster{$key}[2]\t$clustersize\t$cluster{$key}[3]\t$cluster{$key}[4]\t$cluster{$key}[5]\n";
	}
}
