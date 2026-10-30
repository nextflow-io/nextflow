nextflow.enable.types = true

// methylseq `workflows/methylseq/args.nf` defaults rewritten as option maps (see cli() in ../workflows/args.nf)

include { SampleMeta ; TrimgaloreParams ; BismarkParams ; MethyldackelParams } from './types'

def trimgaloreOpts(s: SampleMeta, p: TrimgaloreParams) -> Map<String,?> {
    // Clip presets per protocol: [clip_r1, clip_r2, three_prime_clip_r1, three_prime_clip_r2], 0 = none
    def preset: List<Integer> =
        p.skip_trimming_presets ? [0, 0, 0, 0] :
        p.pbat                  ? [8, 8, 8, 8] :
        p.single_cell           ? [6, 6, 6, 6] :
        p.zymo || p.em_seq      ? [10, 10, 10, 10] :
        p.accel                 ? [10, 15, 10, 10] :
                                  [0, 0, 0, 0]
    return [
        fastqc: true,
        rrbs: p.rrbs,
        nextseq: p.nextseq_trim > 0 ? p.nextseq_trim : null,
        length: p.length_trim ?: null,
        clip_r1: p.clip_r1 > 0 ? p.clip_r1 : preset[0] ?: null,
        clip_r2: s.single_end ? null : p.clip_r2 > 0 ? p.clip_r2 : preset[1] ?: null,
        three_prime_clip_r1: p.three_prime_clip_r1 > 0 ? p.three_prime_clip_r1 : preset[2] ?: null,
        three_prime_clip_r2: s.single_end ? null : p.three_prime_clip_r2 > 0 ? p.three_prime_clip_r2 : preset[3] ?: null
    ]
}

def bismarkAlignOpts(s: SampleMeta, p: BismarkParams) -> Map<String,?> {
    // Combined-index alignment is incompatible with --local_alignment, so gated off there
    def hisat = p.aligner == 'bismark_hisat'
    def non_directional = p.single_cell || p.non_directional || p.zymo
    def use_combined = p.aligner.startsWith('bismark') && p.combined_index && !p.local_alignment
    return [
        hisat2: hisat,
        bowtie2: !hisat,
        'known-splicesite-infile': hisat && p.known_splices ? "<(hisat2_extract_splice_sites.py ${p.known_splices})" : null,
        pbat: p.pbat,
        non_directional: non_directional,
        combined_index: use_combined,
        combined_index_sequential: use_combined && non_directional,
        unmapped: p.unmapped,
        score_min: p.relax_mismatches ? "L,0,-${p.num_mismatches}" : null,
        local: p.local_alignment,
        minins: s.single_end ? null : p.minins,
        maxins: s.single_end ? null : p.maxins ?: (p.em_seq ? 1000 : null)
    ]
}

def bismarkMethylationExtractorOpts(s: SampleMeta, p: BismarkParams) -> Map<String,?> {
    def pe = !s.single_end
    return [
        comprehensive: p.comprehensive,
        cutoff: p.meth_cutoff,
        CX: p.nomeseq,
        ignore: p.ignore_r1 > 0 ? p.ignore_r1 : null,
        ignore_3prime: p.ignore_3prime_r1 > 0 ? p.ignore_3prime_r1 : null,
        no_overlap: pe && p.no_overlap,
        include_overlap: pe && !p.no_overlap,
        ignore_r2: pe && p.ignore_r2 > 0 ? p.ignore_r2 : null,
        ignore_3prime_r2: pe && p.ignore_3prime_r2 > 0 ? p.ignore_3prime_r2 : null
    ]
}

def methyldackelExtractOpts(p: MethyldackelParams) -> Map<String,?> {
    return [
        CHG: p.all_contexts,
        CHH: p.all_contexts,
        mergeContext: p.merge_context,
        ignoreFlags: p.ignore_flags,
        methylKit: p.methyl_kit,
        minDepth: p.min_depth > 0 ? p.min_depth : null
    ]
}
