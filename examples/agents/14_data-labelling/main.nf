nextflow.enable.types = true

// Agentic map-reduce for DATA LABELLING: harmonize messy public sample metadata
// into controlled-vocabulary ontology labels.
//
//   SHARD  (deterministic) : META_FETCH pulls one run's raw metadata from ENA.
//   MAP    (agentic)       : `classify` labels each run's free-text metadata against a
//                            controlled vocabulary (the `metadata-ontology` skill),
//                            emitting one structured label record per sample.
//   REDUCE (agentic)       : `audit` gathers ALL labels, unifies equivalent terms across
//                            the cohort, and splits confident labels from a human-review
//                            queue — emitting a harmonized table.
//
// The deterministic scaffolding (channel = work queue, process fan-out = shard,
// `collect()` = gather) is plain Nextflow; the agents supply only the judgment.
// The `skill` carries the ontology, so labels conform to a fixed vocabulary the model
// would not otherwise use. See README.md.

// One sample's controlled-vocabulary labels, produced by the MAP agent (structured output).
// Each facet is `term [ONTOLOGY:ID]` or the literal `unknown` (see the skill).
record SampleLabels {
    accession:  String    // run accession (from the metadata) — for traceability
    organism:   String    // e.g. "Homo sapiens [NCBITaxon:9606]"
    tissue:     String    // tissue / cell type, e.g. "B lymphocyte [CL:0000236]" or "unknown"
    disease:    String    // e.g. "none [—]" or "unknown"
    assay:      String    // controlled assay term, e.g. "RNA-seq", "ATAC-seq"
    confidence: Double     // 0..1 overall labelling confidence — the review gate keys on this
    rationale:  String    // free-text phrases keyed on; names any facet left `unknown`
}

// The harmonized cohort, produced by the REDUCE agent (structured output).
record HarmonizedCohort {
    tsv:          String  // harmonized table: accession, organism, tissue, disease, assay, confidence
    review_queue: String  // samples/facets below the confidence floor or in conflict — for a human
    summary:      String  // brief: how many labelled confidently, terms unified, conflicts found
}

// SHARD: fetch one run's raw metadata from the ENA "filereport" API. Deterministic,
// no reasoning — runs as a canonical container task on any executor. One task per
// accession (fan-out).
process META_FETCH {
    conda 'conda-forge::curl'

    input:
    accession: String

    output:
    meta: String = stdout()

    script:
        def fields = 'run_accession,sample_accession,sample_title,scientific_name,tax_id,' +
                     'library_strategy,library_source,library_selection,instrument_platform,read_count'
        """
        curl --fail --silent --show-error --location --get \
            --data-urlencode 'accession=${accession}' \
            --data-urlencode 'result=read_run' \
            --data-urlencode 'fields=${fields}' \
            --data-urlencode 'format=json' \
            'https://www.ebi.ac.uk/ena/portal/api/filereport'
        """
}

// MAP: one run's messy metadata -> one controlled-vocabulary label record.
// Structured output + a `skill` (the ontology): the skill loop runs, then a final
// structuring turn encodes the answer into the SampleLabels schema.
agent classify {
    model 'openai/gpt-5-mini'
    instruction '''\
        You are a biomedical metadata curator. Given one sequencing run's raw metadata
        (JSON from the ENA "filereport" API — an array holding one run record), assign
        controlled-vocabulary labels for organism, tissue/cell type, disease, and assay.
        Use the metadata-ontology skill for the allowed terms, ID format, and labelling
        rules. Base every label on evidence in the metadata; when a facet is not
        evidenced, label it `unknown` rather than guessing, and lower the overall
        confidence accordingly. Carry the run accession through for traceability.
        '''.stripIndent()

    // The ontology lives in the skill, not in this prompt — so labels conform to a fixed
    // vocabulary and the same skill can be swapped for an institution's real one.
    skills 'metadata-ontology'

    input:
    meta: String
    output:
    labels: SampleLabels

    prompt:
    """
    Label this sequencing run's metadata:

    ${meta}
    """
}

// REDUCE: all label records -> one harmonized cohort + review queue. Structured, no tools.
agent audit {
    model 'openai/gpt-5-mini'
    instruction '''\
        You reconcile a cohort of per-sample metadata labels into one harmonized table.
          - Unify labels that mean the same thing across samples (same cell line or term
            written differently) so the cohort is internally consistent.
          - Split the cohort: samples whose overall confidence is below 0.6, or that
            disagree with an otherwise-consistent group, go into a human-review queue with
            a one-line reason; the rest are considered accepted.
          - Emit `tsv` as a tab-separated table with the exact header
            `accession\torganism\ttissue\tdisease\tassay\tconfidence` and one row per
            sample, ordered by accession.
          - In `review_queue`, list each flagged sample as `accession: <reason>` (or
            "none" if every sample is confidently labelled).
          - In `summary`, briefly state how many samples were accepted vs queued, which
            terms you unified, and any conflicts.
        '''.stripIndent()

    input:
    cohort: List<SampleLabels>   // fan-in of all labels, deterministically ordered (see toSortedList below)
    output:
    report: HarmonizedCohort

    prompt:
    """
    Reconcile these per-sample labels into a harmonized cohort:

    ${groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(cohort))}
    """
}

workflow {
    // Work queue: real ENA run accessions chosen to exercise the messy cases (edit to
    // point at your own study). Deliberately mixed:
    //   SRR3192396  human  RNA-seq   RICH title: "B-lymphocyte, lymphoblastoid, ...EBV" -> confident
    //   SRR031708   fly    RNA-seq   "S2_DRSC_Untreated-1" -> organism/cell-line from a name
    //   SRR891268   human  ATAC-seq  library_strategy=OTHER but title says ATACseq -> assay from free text
    //   SRR307898   human  RNA-seq   GM12878 again (same line as 891268) -> REDUCE unifies
    //   SRR1039508  human  RNA-seq   "N61311_untreated" -> tissue/disease unevidenced -> review queue
    //   SRR5150592  human  RNA-seq   "Set7KD_rep1" -> unevidenced -> review queue
    def accessions = channel.of(
        'SRR3192396', 'SRR031708', 'SRR891268',
        'SRR307898', 'SRR1039508', 'SRR5150592'
    )

    def meta   = META_FETCH(accessions)   // SHARD  (parallel process fan-out)
    // MAP: even though `classify` declares `skills`, it lowers to a real parallel TaskRun (skills
    // carry no dataflow coupling), fanning out up to executor.cpus with progress rows and resume
    // caching — exactly like a tool-free map agent (cf. samplesheet-builder). Only `tools` agents
    // still take the legacy serial path.
    def labels = classify(meta)

    // Gather all labels into ONE list for the reduce. `toSortedList` (vs `collect()`) fixes the
    // fan-in order by accession, so the output TSV is deterministically ordered. (With `classify`
    // on the task path its output is a task, so the reduce upstream has a stable resume hash and
    // caches too. The expensive map calls are what -resume saves.)
    audit(labels.toSortedList { it.accession }).view { r ->   // REDUCE (gather -> one agent call)
        file("${projectDir}/harmonized.tsv").text = r.tsv
        """\
        Wrote ${projectDir}/harmonized.tsv

        ${r.tsv}
        --- review queue ---
        ${r.review_queue}
        --- summary ---
        ${r.summary}
        """.stripIndent()
    }
}
