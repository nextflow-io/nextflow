nextflow.enable.types = true

record SampleMeta {
    id: String
    single_end: Boolean
    tool_args: Map<String,String>
}
record TrimgaloreParams {
    rrbs: Boolean
    nextseq_trim: Integer
    length_trim: Integer?
    clip_r1: Integer
    clip_r2: Integer
    three_prime_clip_r1: Integer
    three_prime_clip_r2: Integer
    skip_trimming_presets: Boolean
    pbat: Boolean
    single_cell: Boolean
    accel: Boolean
    zymo: Boolean
    em_seq: Boolean
}

record BismarkParams {
    args: Map<String,String>
    aligner: String
    known_splices: Path?
    pbat: Boolean
    single_cell: Boolean
    non_directional: Boolean
    zymo: Boolean
    em_seq: Boolean
    slamseq: Boolean
    combined_index: Boolean
    local_alignment: Boolean
    unmapped: Boolean
    relax_mismatches: Boolean
    num_mismatches: Float
    minins: Integer?
    maxins: Integer?
    comprehensive: Boolean
    meth_cutoff: Integer?
    nomeseq: Boolean
    ignore_r1: Integer
    ignore_3prime_r1: Integer
    no_overlap: Boolean
    ignore_r2: Integer
    ignore_3prime_r2: Integer
}

record MethyldackelParams {
    args: Map<String,String>
    all_contexts: Boolean
    merge_context: Boolean
    ignore_flags: Boolean
    methyl_kit: Boolean
    min_depth: Integer
}

record QualimapParams {
    genome: String?
}

record MethuratorParams {
    methurator_compute_ci: Boolean
    rrbs: Boolean
    methurator_minimum_coverage: String?
    methurator_t_max: Integer?
}

record MultiqcParams {
    args: Map<String,String>
    multiqc_title: String?
}
