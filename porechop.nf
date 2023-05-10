#! /usr/bin/env nextflow

workflow {
   def fastqs = Channel.fromPath (params.pathToFastqs)
   runPoreChop (fastqs)
}

process runPoreChop {

memory '60 GB'

publishDir "$params.porechop_out"

input:
	path fastq
output:
	path "${fastq}_filtered.fastq"
	path "${fastq}_trimmed.fastq"
	
script:
	
"""
seqkit seq ${fastq} -m 100 > ${fastq}_filtered.fastq
wait
porechop -t $params.threads -i ${fastq}_filtered.fastq -o ${fastq}_trimmed.fastq
wait
rm ${fastq}_filtered.fastq
"""

}
