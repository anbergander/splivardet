#! /usr/bin/env nextflow

include { splicevariantdet_wf } from './modules/spliceVariantDetection.nf'


workflow{
   
   def fastqs_trimmed = Channel.fromPath (params.dirToFastqs+'/*.final.fastq')
   

   splicevariantdet_wf (fastqs_trimmed)
}
