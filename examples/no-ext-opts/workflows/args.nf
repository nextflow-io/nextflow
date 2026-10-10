nextflow.enable.types = true

include { Sample ; MiniParams } from './types'

/*
 * Resolve tool args:
 * `<tool>_args` samplesheet column (raw string) > pipeline defaults merged with `--opts.<tool>.<option>`
 */
def toolArgs(tool: String, s: Sample, p: MiniParams, defaults: Map<String,?>) -> String {
    return s.tool_args[tool] ?: cli(defaults + cliOpts(p.opts[tool] ?: [:]))
}

/*
 * Render tool options as CLI args. Boolean -> bare flag or omitted, null -> omitted,
 * any other value -> `<flag> <value>`. The flag is `-k` for a 1-char key, `--key` otherwise,
 * or the key itself if it starts with `-`.
 */
def cli(opts: Map<String,?>) -> String {
    return opts.keySet()
        .collect { k ->
            def v = opts[k]
            def flag = k.startsWith('-') ? k : (k.length() == 1 ? "-${k}" : "--${k}")
            v instanceof Boolean ? (v ? flag : '') : v != null ? "${flag} ${v}" : ''
        }
        .findAll { a -> a != '' }
        .join(' ')
}

// CLI values arrive as strings: 'true'/'false' mean flag on/off
def cliOpts(opts: Map<String,?>) -> Map<String,?> {
    return opts.keySet().inject([:]) { acc, k ->
        def v = "${opts[k]}"
        acc + [(k): opts[k]] + (v == 'true' ? [(k): true] : v == 'false' ? [(k): false] : [:])
    }
}

def trimgaloreOpts(s: Sample, p: MiniParams) -> Map<String,?> {
    return [
        fastqc: true,
        rrbs: p.rrbs,
        clip_r1: s.clip_r1 > 0 ? s.clip_r1 : (p.skip_trimming_presets ? null : 10),
        clip_r2: !s.single_end && p.clip_r2 > 0 ? p.clip_r2 : null
    ]
}

def alignOpts(s: Sample, p: MiniParams) -> Map<String,?> {
    return [
        bowtie2: true,
        maxins: s.single_end ? null : p.maxins
    ]
}
