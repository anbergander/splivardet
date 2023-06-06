#! /usr/bin/env nextflow

workflow splicevariantdet_wf{
take:
	fastqs_trimmed

main:
	runDesaltMapping (fastqs_trimmed)
	runGenerateBamBed (runDesaltMapping.out)
	runIsoformDetectionFlair (runGenerateBamBed.out)

}



process runDesaltMapping

{

memory '60 GB'

input:
	path fastq


output:
	path "${fastq}.sam"
	path "${fastq}"

script:

"""

deSALT aln -s 2 -t $params.threads -l 14 -x ont1d $params.desaltIndex -f ${fastq}.temp $fastq -o ${fastq}.sam

 
exit 0
"""

}

process runGenerateBamBed
{

input:
	path sam
	path fastq
	
output:
	path "${sam}.bam"
	path "${sam}.bam.bai"
	path "${sam}.bed"
	path "${fastq}"
	
script:


"""
samtools view -bS $sam|samtools sort -@ $params.threads -o ${sam}.bam

samtools index ${sam}.bam -@ $params.threads ${sam}.bam.bai 

bedtools bamtobed -bed12 -i ${sam}.bam > ${sam}.bed
 
exit 0
"""

}

process runIsoformDetectionFlair
{
publishDir "$params.out",mode:"move"
input:
	path bam
	path bambai
	path bed
	path fastq

	
output:
	path "${out}_all_corrected.bed"
	path "${out}_all_inconsistent.bed", optional: true
	path "${out}_cannot_verify.bed", optional: true
	path "${out}.isoforms.bed"
	path "${out}.isoforms.fa"
	path "${out_fq}.final.fastq"
	
script:

out = bed.getSimpleName()
out_fq = fastq.getSimpleName()

"""
flair correct -q $bed -f $params.annotationsGtf -g $params.genomeFasta -o ${out}

flair collapse -g $params.genomeFasta -q ${out}_all_corrected.bed -r $fastq -o ${out}

exit 0
"""

}

