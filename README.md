# PTC_predictor

Perl script for calculating location of premature termination codon for coding frameshift mutations. Tested on Linux and WSL.

## Prerequisites
> 1. Linux, Mac OS, or WSL
> 2. Perl
> 3. Ensembl Perl API (including Bioperl)
> 4. GNU parallel (optional)

## Run on single position
Simply clone git directory and provide sample input. No need to install or make files.
````shell
git clone https://github.com/ITSpaanderman/PTC_predictor.git
cd ./PTC_predictor
perl PTC_predictor [ROWNMUBER] [FIELD]
# FIELD conisist of mutation locus, muttype, mutsubtype, mutchange
# example: chr1:123456-INS-Frameshift-1
````
Locus: chr[n]:[position] <br>
Muttypes: INS, DEL, SNP <br>
Mutsubtypes: Frameshift, In_Frame, Nonsense, SNP, None, Undef <br>
Mutchange: [n] <br>

## Run in parallel
Clone git and provide input file with rownumbers and field information conisiting of locus, muttype, mutsubtpye, mutchange (see above)
````shell
git clone https://github.com/ITSpaanderman/PTC_predictor.git
cd ./PTC_predictor
cat [input_file.csv] | parallel --colsep ',' PTC_predictor.pl {1} {2} > [output_file.tsv]
````

## Test
PTC_test_in.tsv input file and corresponding PTC_test_out.tsv output file are located in Test folder.
