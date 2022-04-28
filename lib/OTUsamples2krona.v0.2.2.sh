#!/bin/bash
NAME='OTUsamples2krona'
VER='v0.2.2' # version for Github
AUTHOR='Bruno Gomez-Gil' # Laboratorio de Genomica Microbiana, CIAD. https://github.com/GenomicaMicrob
REV='May 21, 2018'
DEP='KronaTools v2.5'
LINK='https://github.com/GenomicaMicrob/OTUsamples2krona'

display_version(){
echo "
____________________________________________________________________
Script name:    $NAME
Version:        $VER
Last revisited: $REV
Dependencies:   $DEP
Author:         $AUTHOR
More info at:   $LINK
____________________________________________________________________
"
} # If -v is typed, displays version
if [ "$1" == "-v" ]
	then
		display_version
		exit 1
fi
display_help(){
echo -e "
_____________________________ $NAME $VER _________________________________________________________________
A simple script to generate Krona pie charts from all samples OTUs. It also produces a file with the OTU,
  number of hits, and percentage.
\e[1mUSAGE\e[0m: $NAME file.tsv
The tsv file has to have the following format:
   first line has the sample name (delimited by tabs) and taxa structure (delimited by semicolons)
   from the second line, the OTU values (delimited by tabs) and the taxonomy delimited by semicolons.
   
 Example:
 sample1 sample2 sampleN domain;phylum;class;order;family;genus;species
   2       0       90    Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Marinococcus;Unclassified_s
   12      300     0     Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus_iners
mg_classifier produces already the file samples-tax.tsv with this format.
_______________________________________________________________________________________________________________________
"
} # -h is typed, displays help
if [ "$1" == "-h" ]
	then
		display_help
		exit 1
fi
display_usage(){
echo -e "
__________________ $NAME $VER _________________________
\e[1mERROR\e[0m: missing file
\e[1mUSAGE\e[0m: $NAME samples-tax.tsv
For help, type: $NAME -h
______________________________________________________________________
"
} # less than one arguments supplied, displays usage 
if [  $# -le 0 ] 
	then 
		display_usage
		exit 1
fi
# ------- script starts -------
NCOL=$(awk 'NR==1{print NF-1}' $1) # variable to get the number of columns, except the last (taxonomy)
mkdir -p krona
awk '{print $NF}' $1 > krona/tax.txt # creates a file with only the taxonomy for later pasting
for (( i=1; i <= $NCOL; i++ )); do awk '{print $'$i'}' $1 > krona/sample_0${i}.dat; done # makes temporary files with the info of each sample (column)
# create individual files with OTU values and taxonomy

for f in krona/*.dat
do
	SAMPLE=$(head -1 $f)
	paste $f krona/tax.txt > krona/$SAMPLE.tmp
	# Generating file with the otus per sample
	paste krona/tax.txt $f > krona/$SAMPLE.id
	grep -w -v "0" krona/$SAMPLE.id | sed 1d | awk '{sums[$1] += $2} END { for (i in sums) printf("%s %s\n", i, sums[i])}' | sed 's/ /\t/' | sort -n -r -k2,2 > krona/$SAMPLE.txt
	cut -f2 krona/$SAMPLE.txt | awk '{array[NR] = $0; sum+= $0 } END {for (x = 1; x <= NR; x++) printf "%2.2f\n", (100 * array[x])/sum }' > krona/$SAMPLE.percent
	paste krona/$SAMPLE.txt krona/$SAMPLE.percent > krona/$SAMPLE.tsv
	sed -i '1i domain;phylum;class;order;family;genus;species\thits\tpercent' krona/$SAMPLE.tsv
done # formatting files for krona

for f in krona/*.tmp
do
	sed -i '1d' $f
	grep -w -v "0" $f > $f.tax
done # removes first line (header) and lines with zero as value

rename 's/tmp.//' krona/*.tax # removes tmp. from all .tax file names

for f in krona/*.tax
do
	sort -k2,2 $f -o $f # sorts based on taxonomy
	sed -i 's/;/\t/g' $f
	#SAMPLE=$(basename $f .tax) # gets the filename without the extension
	SAMPLE=$(basename -s .tax $f ) # gets the filename without the extension
	ktImportText $f -n $SAMPLE -o $f.html # generates Krona charts
done # formatting and charting

# clean up
rename -f 's/tax.//' krona/*.html
rm -f krona/*.dat krona/*.tmp krona/*.txt krona/*.tax krona/*.id krona/*.percent
#[edited by]echo -e "\e[33m$NCOL Krona charts have been created in krona/\e[0m"
echo -e "Krona charts have been created in krona"
echo "Adios"
echo
# This is the end.
