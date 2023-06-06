#! /usr/bin/env nextflow

include { porechop_wf } from './modules/porechop.nf'


workflow{
   def fq_gz_ch = Channel.fromPath (params.dirToFastqs+'*.fastq.gz')
   

   porechop_wf (fq_gz_ch)
}
