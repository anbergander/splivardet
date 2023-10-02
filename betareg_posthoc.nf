#! /usr/bin/env nextflow

workflow{

def fastqs = Channel.fromPath (params.dirToFastqs+'*.final.fastq')
def quantfiles = Channel.fromPath (params.out+'/*.quantify*')

countReads(fastqs)
runStatandPlot(quantfiles)

}

process cleanup{


"""
[ ! -e ${params.out}read_counts.csv ] || rm ${params.out}read_counts.csv

[ ! -e ${params.out}items.txt ] || rm ${params.out}items.*.txt

"""

}


process runStatandPlot{

publishDir "$params.out",mode:"move"

input:
	path quantfile
output:
	path "*.pdf"
	path "${b}.candidates.txt"
script:

b = quantfile.getSimpleName()

"""
python3 $workflow.projectDir/preprocess_stat.py $quantfile ${params.out}read_counts.csv
Rscript $workflow.projectDir/betareg.R
python3 $workflow.projectDir/fdr_correction.py
python3 $workflow.projectDir/generate_plots.py
mv candidates.txt ${b}.candidates.txt
rm *pickle

"""

}

process countReads {

publishDir "$params.out",mode:"move"

input:
	path fastq
	
script:

b = fastq.getSimpleName()

"""
[ -e  ${params.out}/read_counts.csv] || touch ${params.out}/read_counts.csv
echo ${b},\$(echo \$(cat $fastq|wc -l)/4|bc) >> ${params.out}/read_counts.csv
"""


}


