#! /usr/bin/env nextflow

workflow porechop_wf{
take:
   fastqs_ch
   
main:
   runPoreChop (fastqs_ch)
  
}

process runPoreChop {

memory '40 GB'

publishDir "$params.out",mode:"move"

input:
	path fastq
output:
	path "${b}.final.fastq", optional:true
	path "${b}.trimmed.fastq.gz", optional:true
	
	

script:


b = fastq.getSimpleName()

	
"""
porechop -t $params.threads -i ${fastq} -o ${b}.trimmed.fastq.gz
seqkit seq ${b}.trimmed.fastq.gz -m 100 > ${b}.final.fastq
rm ${b}.trimmed.fastq.gz

"""

}
