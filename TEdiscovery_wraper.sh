#!/bin/sh

# This script calls scripts in order to identify new elements

TEdiscoveryPath="../TEdiscovery" # path to scripts
gname=$1 # nickname of genome
consensi=$2 # repeatmodeler consensi
genome=$3 # genome fasta file
dbgenome=$4 # blast db for genome

# compare the repeatmodeler output to known elements
consensi2outputfile=$gname-consensi2.fa
perl $TEdiscoveryPath/search_potential_elements.pl -t $TEdiscoveryPath/repbase_classII.fa -p $consensi -g $genome -d $dbgenome -o $consensi2outputfile 

# run repeatmasker
~/bin/RepeatMasker/RepeatMasker -lib $consensi2outputfile -dir . $genome -pa 12 &
