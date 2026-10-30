nextflow.enable.types = true

include { Sample ; MiniParams } from './types'

// Resolve tool args: samplesheet column > --args.<tool> > pipeline default
def toolArgs(tool: String, s: Sample, p: MiniParams, fallback: String) -> String {
    return s.tool_args[tool] ?: p.args[tool] ?: fallback
}

def trimgaloreArgs(s: Sample, p: MiniParams) -> String {
    return [
        '--fastqc',
        p.rrbs ? '--rrbs' : '',
        s.clip_r1 > 0 ? "--clip_r1 ${s.clip_r1}" : (p.skip_trimming_presets ? '' : '--clip_r1 10'),
        s.single_end ? '' : (p.clip_r2 > 0 ? "--clip_r2 ${p.clip_r2}" : '')
    ].findAll { a -> a != '' }.join(' ')
}

def alignArgs(s: Sample, p: MiniParams) -> String {
    return [
        '--bowtie2',
        !s.single_end && p.maxins ? "--maxins ${p.maxins}" : ''
    ].findAll { a -> a != '' }.join(' ')
}
