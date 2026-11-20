nextflow.enable.types = true

// Asserts the option-map defaults render the same args as the string versions, for every combination of the params they read

include { trimgaloreArgs ; bismarkAlignArgs ; bismarkMethylationExtractorArgs } from './strings'
include { trimgaloreOpts ; bismarkAlignOpts ; bismarkMethylationExtractorOpts } from './opts'
include { cli } from '../workflows/args'

def compare(tool: String, str: String, opts: String) -> Integer {
    def a = str.tokenize(' ').join(' ')
    def b = opts.tokenize(' ').join(' ')
    if( a == b )
        return 0
    println "${tool}\n  str:  ${a}\n  opts: ${b}"
    return 1
}

workflow {
    def B = [true, false]
    def trim = [B, B, B, B, B, B, B, B, [0, 3], [null, 20]].combinations().collect { c ->
        def s = record(id: 's', single_end: c[0], tool_args: [:])
        def p = record(rrbs: c[1], pbat: c[2], single_cell: c[3], accel: c[4], zymo: c[5], em_seq: c[6], skip_trimming_presets: c[7],
            nextseq_trim: c[8], clip_r1: c[8], clip_r2: c[8], three_prime_clip_r1: 0, three_prime_clip_r2: c[8], length_trim: c[9])
        compare('trimgalore', trimgaloreArgs(s, p), cli(trimgaloreOpts(s, p)))
    }
    def bismark = [B, ['bismark', 'bismark_hisat', 'bwameth'], [null, file('splices.txt')], B, B, B, B, B, B, B, B, [null, 50], [null, 500], B].combinations().collect { c ->
        def s = record(id: 's', single_end: c[0], tool_args: [:])
        def p = record(aligner: c[1], known_splices: c[2], pbat: c[3], single_cell: c[4], non_directional: c[5], zymo: c[6],
            combined_index: c[7], local_alignment: c[8], unmapped: c[9], relax_mismatches: c[10], num_mismatches: 0.6, minins: c[11], maxins: c[12], em_seq: c[13],
            comprehensive: c[3], meth_cutoff: c[11], nomeseq: c[4], ignore_r1: c[5] ? 2 : 0, ignore_3prime_r1: c[6] ? 2 : 0,
            no_overlap: c[7], ignore_r2: c[8] ? 2 : 0, ignore_3prime_r2: c[9] ? 2 : 0)
        compare('bismark_align', bismarkAlignArgs(s, p), cli(bismarkAlignOpts(s, p)))
            + compare('methylation_extractor', bismarkMethylationExtractorArgs(s, p), cli(bismarkMethylationExtractorOpts(s, p)))
    }
    println "trimgalore: ${trim.size()} cases, ${trim.sum()} mismatches"
    println "bismark align + methylation extractor: ${bismark.size() * 2} cases, ${bismark.sum()} mismatches"
}
