#! /usr/bin/env nextflow

workflow quantify_wf{
take:
isoformFastas

   
main:
generateDiffManifest(params.flairManifest)
runFlairQuantify(isoformFastas, params.flairManifest)



}



process generateDiffManifest {

input:
	path manifest
output:
	path 'diff.manifest.tab'
	
shell:

"""
awk -F '\\t' '{if (\$2 == "$params.condition_group"||\$2 == "$params.control_group")print \$0}' $manifest > diff.manifest.tab
"""


}




process runFlairQuantify {

publishDir "${params.out}", mode:'move'

input:
	path isoformFasta
	path diffManifest

output:
	path "${b}.quantify.counts.tsv"
maxForks=1

script:


b = isoformFasta.getSimpleName()


	
"""
flair quantify -r $diffManifest -i $isoformFasta --quality 4 --temp_dir ${params.out}temp --output ${b}.quantify
rm -r ${params.out}temp

"""

}
