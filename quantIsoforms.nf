#! /usr/bin/env nextflow

include { quantify_wf; generateDiffManifest } from './modules/quantifySpliceVariants.nf'


workflow{
   def isoform_ch = Channel.fromPath (params.dirToFastqs+'*.fa')
   
   
   quantify_wf (isoform_ch)
}
