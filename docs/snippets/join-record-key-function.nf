reads = channel.of( record(meta: [id: 'X'], fastq: 'X.fq'), record(meta: [id: 'Y'], fastq: 'Y.fq') )
bams  = channel.of( record(meta: [id: 'Y', strandedness: 'reverse'], bam: 'Y.bam'), record(meta: [id: 'X', strandedness: 'forward'], bam: 'X.bam') )

reads.join(bams, by: { r -> r.meta.id }).view()
