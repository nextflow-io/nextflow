nextflow.preview.types = true

workflow {
    channel.of(
            tuple('chr1', ['/path/to/region1_chr1.vcf', '/path/to/region2_chr1.vcf']),
            tuple('chr2', ['/path/to/region1_chr2.vcf', '/path/to/region2_chr2.vcf', '/path/to/region3_chr2.vcf']),
        )
        .flatMap { chr, vcfs ->
            vcfs.collect { vcf ->
                tuple(chr, vcfs.size(), vcf)    // preserve group size
            }
        }
        .view { v -> "scattered: ${v}" }
        .groupBy()
        .map { key, values -> tuple(key, values.toSorted()) }
        .view { v -> "gathered: ${v}" }
}