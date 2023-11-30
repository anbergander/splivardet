#! /usr/bin/env nextflow

workflow{

def gtfs = Channel.fromPath (params.out+'*.isoforms.gtf')

runFormatGtf(gtfs)

}

process runFormatGtf{

publishDir "$params.out",mode:"move"

input:
	path gtf

output:
	path "${b}.sign.gtf"
script:

b = gtf.getSimpleName()

"""
python3 $workflow.projectDir/format_isoform_gtf.py $gtf ${params.out}${b}.candidates.txt


"""

}
