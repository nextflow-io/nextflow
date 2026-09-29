nextflow.enable.types = true

// Verbatim from methylseq workflows/methylseq/args.nf (string-building versions)

include { SampleMeta ; TrimgaloreParams ; BismarkParams } from './types'

def trimgaloreArgs(s: SampleMeta, p: TrimgaloreParams) -> String {
    return [
        // Static args
        '--fastqc',

        // Special flags
        p.rrbs ? '--rrbs' : '',
        p.nextseq_trim > 0 ? "--nextseq ${p.nextseq_trim}" : '',
        p.length_trim ? "--length ${p.length_trim}" : '',

        // Trimming - R1
        p.clip_r1 > 0 ? "--clip_r1 ${p.clip_r1}" : (
            p.skip_trimming_presets ? '' : (
                p.pbat ? "--clip_r1 8" : (
                    p.single_cell ? "--clip_r1 6" : (
                        (p.accel || p.zymo || p.em_seq) ? "--clip_r1 10" : ''
                    )
                )
            )
        ),

        // Trimming - R2
        s.single_end ? '' : (
            p.clip_r2 > 0 ? "--clip_r2 ${p.clip_r2}" : (
                p.skip_trimming_presets ? '' : (
                    p.pbat ? "--clip_r2 8" : (
                        p.single_cell ? "--clip_r2 6" : (
                            (p.zymo || p.em_seq) ? "--clip_r2 10" : (
                                p.accel ? "--clip_r2 15" : ''
                            )
                        )
                    )
                )
            )
        ),

        // Trimming - 3' R1
        p.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${p.three_prime_clip_r1}" : (
            p.skip_trimming_presets ? '' : (
                p.pbat ? "--three_prime_clip_r1 8" : (
                    p.single_cell ? "--three_prime_clip_r1 6" : (
                        (p.accel || p.zymo || p.em_seq) ? "--three_prime_clip_r1 10" : ''
                    )
                )
            )
        ),

        // Trimming - 3' R2
        s.single_end ? '' : (
            p.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${p.three_prime_clip_r2}" : (
                p.skip_trimming_presets ? '' : (
                    p.pbat ? "--three_prime_clip_r2 8" : (
                        p.single_cell ? "--three_prime_clip_r2 6" : (
                            (p.accel || p.zymo || p.em_seq) ? "--three_prime_clip_r2 10" : ''
                        )
                    )
                )
            )
        ),
    ].join(' ').trim()
}

def bismarkAlignArgs(s: SampleMeta, p: BismarkParams) -> String {
    // Combined-index alignment is incompatible with --local_alignment, so gated off there
    def non_directional = p.single_cell || p.non_directional || p.zymo
    def use_combined = p.aligner.startsWith('bismark') && p.combined_index && !p.local_alignment
    return [
        (p.aligner == 'bismark_hisat') ? ' --hisat2' : ' --bowtie2',
        (p.aligner == 'bismark_hisat' && p.known_splices) ? " --known-splicesite-infile <(hisat2_extract_splice_sites.py ${p.known_splices})" : '',
        p.pbat ? ' --pbat' : '',
        non_directional ? ' --non_directional' : '',
        use_combined ? ' --combined_index' : '',
        (use_combined && non_directional) ? ' --combined_index_sequential' : '',
        p.unmapped ? ' --unmapped' : '',
        p.relax_mismatches ? " --score_min L,0,-${p.num_mismatches}" : '',
        p.local_alignment ? " --local" : '',
        !s.single_end && p.minins ? " --minins ${p.minins}" : '',
        s.single_end ? '' : (
            p.maxins ? " --maxins ${p.maxins}" : (
                p.em_seq ? " --maxins 1000" : ''
            )
        )
    ].join(' ').trim()
}

def bismarkMethylationExtractorArgs(s: SampleMeta, p: BismarkParams) -> String {
    return [
        p.comprehensive   ? ' --comprehensive' : '',
        p.meth_cutoff     ? " --cutoff ${p.meth_cutoff}" : '',
        p.nomeseq         ? '--CX' : '',
        p.ignore_r1 > 0   ? "--ignore ${p.ignore_r1}" : '',
        p.ignore_3prime_r1 > 0   ? "--ignore_3prime ${p.ignore_3prime_r1}" : '',
        s.single_end ? '' : (p.no_overlap           ? ' --no_overlap'                         : '--include_overlap'),
        s.single_end ? '' : (p.ignore_r2        > 0 ? "--ignore_r2 ${p.ignore_r2}"       : ""),
        s.single_end ? '' : (p.ignore_3prime_r2 > 0 ? "--ignore_3prime_r2 ${p.ignore_3prime_r2}": "")
    ].join(' ').trim()
}

