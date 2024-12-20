params.outdir = 'publishDir'

process gen_data {
	publishDir "${params.outdir}", mode: 'copy'
        input:
	val(i)
        output:
        path 'chunk_*'

    	script:
    	"""
   	dd if=/dev/urandom of=chunk_$i count=1 bs=10M
    	"""
}

workflow {
	Channel.of(1) | gen_data 
}
