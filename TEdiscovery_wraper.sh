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
~/bin/RepeatMasker/RepeatMasker -lib $consensi2outputfile $genome -pa 12 &

# convert repeatmasker to bed and filter
rmoutput=$genome.out
rm2gffout=$gname-ele.bed
perl $TEdiscoveryPath/rm2gff3.pl -r $rmoutput -c $consensi2outputfile -o $rm2gffout

# find the TIR locations
perl $TEdiscoveryPath/id-TA-tirs-genome.pl $genome > $gname-tir-TA.bed
perl $TEdiscoveryPath/id-tirs-genome.pl -g $genome -t 4 -i 11 -o $gname-tir-4TSD.bed
perl $TEdiscoveryPath/id-tirs-genome.pl -g $genome -t 8 -i 11 -o $gname-tir-8TSD.bed

# produce reports
perl $TEdiscoveryPath/report-tir-ele-overlap.pl -r $rm2gffout -t $gname-tir-TA.bed -g $genome -d $TEdiscoveryPath/repbase_classII.fa -o $gname-tirele-TA-report.xls
perl $TEdiscoveryPath/report-tir-ele-overlap.pl -r $rm2gffout -t $gname-tir-4TSD.bed -g $genome -d $TEdiscoveryPath/repbase_classII.fa -o $gname-tirele-4TSD-report.xls 
perl $TEdiscoveryPath/report-tir-ele-overlap.pl -r $rm2gffout -t $gname-tir-8TSD.bed -g $genome -d $TEdiscoveryPath/repbase_classII.fa -o $gname-tirele-8TSD-report.xls

exit
