chr_frequency = ["chr1": 2, "chr2": 3]

Channel.of(
        ['region1', 'chr1', '/path/to/region1_chr1.vcf'],
        ['region2', 'chr1', '/path/to/region2_chr1.vcf'],
        ['region1', 'chr2', '/path/to/region1_chr2.vcf'],
        ['region2', 'chr2', '/path/to/region2_chr2.vcf'],
        ['region3', 'chr2', '/path/to/region3_chr2.vcf']
    )
    .map { region, chr, vcf -> tuple( groupKey(chr, chr_frequency[chr]), vcf ) }
    .groupTuple()
    .view()