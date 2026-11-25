nextflow.enable.types = true

// One row of the samplesheet: columns are fields, no meta map
record Sample {
    id: String
    single_end: Boolean
    reads: List<Path>
    clip_r1: Integer
    tool_args: Map<String,String>
}

// Subset of pipeline params used by MINI
record MiniParams {
    rrbs: Boolean
    clip_r2: Integer
    skip_trimming_presets: Boolean
    skip_dedup: Boolean
    maxins: Integer?
    multiqc_title: String?
    args: Map<String,String>
}
