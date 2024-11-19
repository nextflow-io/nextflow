Channel.of(
        ['chr1', ['/path/to/region1_chr1.vcf', '/path/to/region2_chr1.vcf']],
        ['chr2', ['/path/to/region1_chr2.vcf', '/path/to/region2_chr2.vcf', '/path/to/region3_chr2.vcf']],
    )
    .flatMap { chr, vcfs ->
        vcfs.collect { vcf ->
            tuple(groupKey(chr, vcfs.size()), vcf)              // preserve group size with key
        }
    }
    .view { v -> "scattered: ${v}" }
    .groupTuple()
    .map { key, vcfs -> tuple(key.getGroupTarget(), vcfs) }     // unwrap group key
    .view { v -> "gathered: ${v}" }