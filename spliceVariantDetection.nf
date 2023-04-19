#! /usr/bin/env nextflow

workflow {
   def fastqs = Channel.fromPath (params.pathToFastqs)
   runDesaltMapping (fastqs)
   runGenerateBamBed (runDesaltMapping.out)
   runIsoformDetectionFlair (fastqs, runGenerateBamBed.out)
}



process runDesaltMapping

{

memory '60 GB'

input:
	path fastq


output:
	path "${fastq}.sam"

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
	
output:
	path "${sam}.bam"
	path "${sam}.bam.bai"
	path "${sam}.bed"
	
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

input:
	path fastq
	path bam
	path bambai
	path bed
	
output:
	path "${bed}_all_corrected.bed"
	path "${bed}_all_inconsistent.bed", optional: true
	path "${bed}_cannot_verify.bed", optional: true
	path "${bed}.isoforms.bed"
	path "${bed}.isoforms.fa"
	
script:


"""
	flair correct -q $bed -f $params.annotationsGtf -g $params.genomeFasta -o ${bed}
	flair collapse -g $params.genomeFasta -q ${bed}_all_corrected.bed -r $fastq -o ${bed}

exit 0
"""

}

