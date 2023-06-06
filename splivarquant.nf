#! /usr/bin/env nextflow

include { quantify_wf } from './modules/quantifySpliceVariants.nf'


workflow{
   def isoform_fastas = Channel.fromPath (params.dirToFastqs+'*isoforms.fa')
   

   quantify_wf (isoform_fastas)
}
