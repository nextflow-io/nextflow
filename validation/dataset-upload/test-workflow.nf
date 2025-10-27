#!/usr/bin/env nextflow

/*
 * Test workflow for dataset upload functionality
 *
 * This workflow creates a simple CSV output and tests
 * the automatic upload to Seqera Platform datasets.
 */

workflow {
    // Create test sample data
    Channel.of(
        [id: 'sample1', fastq: 'sample1_R1.fq.gz', fastq2: 'sample1_R2.fq.gz'],
        [id: 'sample2', fastq: 'sample2_R1.fq.gz', fastq2: 'sample2_R2.fq.gz'],
        [id: 'sample3', fastq: 'sample3_R1.fq.gz', fastq2: 'sample3_R2.fq.gz']
    )
    | map { meta ->
        // Simulate processing results
        def result = "${meta.id},${meta.fastq},${meta.fastq2},result_${meta.id}.bam"
        result
    }
    | collectFile(
        name: 'results.csv',
        newLine: true,
        storeDir: 'results',
        seed: 'sample_id,fastq_r1,fastq_r2,bam_file'
    )
    | set { results_ch }

    // Publish channel for output block
    publish:
    analysis_results = results_ch
}

// Define workflow outputs for dataset upload
output {
    analysis_results {
        path '.'
        index {
            path 'results.csv'
            header true
        }
    }
}
