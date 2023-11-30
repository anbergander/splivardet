#! /usr/bin/env nextflow

workflow{

def quantfiles = Channel.fromPath (params.out+'*.quantify*')
def isoformfiles = Channel.fromPath (params.out+'*.isoform*')


runStatandPlot(quantfiles)

}



process runStatandPlot{

publishDir "$params.out",mode:"move"

input:
	path quantfile

script:

b = quantfile.getSimpleName()

"""
python3 $workflow.projectDir/preprocess_stat.py ${params.out}${b}.isoforms.gtf ${params.out}read_counts.csv $quantfile $projectDir
python3 $workflow.projectDir/reformat_swan_in.py ${params.out}${b}.isoforms.gtf $projectDir ${params.indexDir} ${params.out}${b}.candidates.txt
python3 $workflow.projectDir/isoform_report.py ${params.out}$b ${params.indexDir} ${params.flairManifest}
mv candidates.txt ${params.out}$b/

"""

}

process countReads {

publishDir "$params.out",mode:"move"

input:
	path fastq
	
script:

b = fastq.getSimpleName()

"""
[ -e  ${params.out}read_counts.csv] || touch ${params.out}read_counts.csv
echo ${b},\$(echo \$(cat $fastq|wc -l)/4|bc) >> ${params.out}read_counts.csv
"""


}


