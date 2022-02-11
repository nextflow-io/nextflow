#!/usr/bin/env nextflow

params.samples_sheet = "$baseDir/data/scratch_csv_bom.csv"

samples_ch = Channel.fromPath(params.samples_sheet)
	 .splitCsv(header:true, strip:true, charset:'UTF-8')
	 .map{ row -> tuple(row.sampleId, row.read1, row.read2) }

process printInput {
	input:
	set sampleId, read1, read2 from samples_ch

	
	output:
	stdout result
    
	script:
	"""
	echo $sampleId $read1 $read2
    """
}

result.view { it }
