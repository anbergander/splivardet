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
	path "${b}.candidates.txt"
script:

b = quantfile.getSimpleName()

"""
python3 $workflow.projectDir/preprocess_stat.py ${params.out}${b}.isoforms.gtf ${params.out}read_counts.csv $quantfile $projectDir
Rscript $workflow.projectDir/betareg.R
python3 $workflow.projectDir/fdr_correction.py $projectDir
mkdir -p ${params.out}${b}
mv posthoc_result.tsv ${params.out}${b}/${b}.posthoc_result.tsv
python3 $workflow.projectDir/generate_plots.py $projectDir
mv *violin.pdf ${params.out}${b}/
mv candidates.txt ${b}.candidates.txt 


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


