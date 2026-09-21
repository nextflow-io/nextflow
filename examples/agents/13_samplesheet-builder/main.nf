nextflow.enable.types = true

// Agentic map-reduce: build a valid FASTQ samplesheet from messy public metadata.
//
//   SHARD (deterministic) : ENA_FETCH pulls raw run metadata, one task per accession.
//   MAP   (agentic)       : `normalize` turns each run's messy metadata into one clean
//                           samplesheet candidate (picks R1/R2, handles single/paired,
//                           derives a safe sample name).
//   REDUCE (agentic)      : `consolidate` gathers ALL candidates and reconciles them into
//                           one samplesheet — merging technical replicates, guaranteeing
//                           unique sample names, flagging conflicts.
//
// The deterministic scaffolding (channel = work queue, process fan-out = shard,
// `collect()` = gather) is plain Nextflow; the agents supply only the judgment.
// See README.md.

// One normalized samplesheet row, produced by the MAP agent (structured output).
record Candidate {
    sample:    String     // clean, filesystem-safe sample name derived from the metadata
    run:       String     // run accession (traceability)
    biosample: String     // biological-sample accession (so REDUCE can merge replicates)
    fastq_1:   String     // full ftp:// URL of R1 (or the single-end reads)
    fastq_2:   String     // full ftp:// URL of R2, or "" for single-end
    layout:    String     // "single" | "paired" — the agent's determination
    note:      String     // any non-obvious choice the agent made (e.g. dropped an orphan file)
}

// The final samplesheet, produced by the REDUCE agent (structured output).
record Samplesheet {
    csv:   String         // full samplesheet CSV text (header + one row per run)
    notes: String         // human-readable summary of merges / drops / conflicts
}

// SHARD: fetch one run's raw metadata from the ENA "filereport" API. Deterministic,
// no reasoning — runs locally (exec), no container. One task per accession (fan-out).
process ENA_FETCH {
    input:
    accession: String
    output:
    meta: String
    exec:
        def fields = 'run_accession,sample_accession,sample_title,library_layout,' +
                     'library_strategy,library_name,scientific_name,fastq_ftp,read_count'
        def url = "https://www.ebi.ac.uk/ena/portal/api/filereport" +
                  "?accession=${accession}&result=read_run&fields=${fields}&format=json"
        meta = new URL(url).getText('UTF-8')
}

// MAP: one run's messy ENA metadata -> one clean candidate. Structured output, no tools.
agent normalize {
    model 'openai/gpt-5-mini'
    instruction '''\
        You convert ONE sequencing run's raw ENA metadata into a single normalized
        FASTQ samplesheet entry. The input is JSON from the ENA "filereport" API
        (an array holding one run record). Work carefully:
          - Determine the library layout (single- or paired-end) from the metadata.
          - The fastq_ftp field lists download URLs joined by ";". For paired-end runs
            ENA sometimes ALSO lists an extra unpaired/orphan file next to the _1/_2
            mates — use only the two mate files (ending _1 and _2) as fastq_1 and
            fastq_2 and ignore the orphan. For single-end runs use the single file as
            fastq_1 and leave fastq_2 empty.
          - Return every FASTQ as a complete URL with an ftp:// scheme.
          - Derive a short, human-readable, filesystem-safe sample name (letters,
            digits, underscores) from the sample title/description. Name runs of the
            same biological sample consistently so they can be grouped later.
          - Carry the run accession and the biological-sample accession through so
            downstream consolidation can merge technical replicates.
          - If you drop a file or make a non-obvious choice, note it briefly.
        Base every value on the metadata; never invent a URL or an accession.
        '''.stripIndent()

    input:
    meta: String
    output:
    candidate: Candidate

    prompt:
    """
    Normalize this ENA run record into one samplesheet candidate:

    ${meta}
    """
}

// REDUCE: all candidates -> one reconciled samplesheet. Structured output, no tools.
agent consolidate {
    model 'openai/gpt-5-mini'
    instruction '''\
        You assemble a set of normalized run entries into ONE valid nf-core-style
        FASTQ samplesheet. Produce CSV with the exact header `sample,fastq_1,fastq_2`
        and one row per run:
          - Group runs belonging to the same biological sample (same biosample
            accession, or clearly the same sample) under ONE consistent `sample`
            name, so nf-core merges them as technical replicates. Give every distinct
            biological sample a unique, filesystem-safe name; disambiguate collisions.
          - Preserve fastq_1 and fastq_2 exactly as provided (fastq_2 empty for
            single-end runs).
          - Order rows stably, grouped by sample.
        Put the complete CSV text in `csv`. In `notes`, briefly summarize your
        decisions: which runs you merged, files dropped upstream, and any conflicts
        (e.g. a sample mixing single- and paired-end runs).
        '''.stripIndent()

    input:
    candidates: List<Candidate>
    output:
    sheet: Samplesheet

    prompt:
    """
    Assemble the samplesheet from these normalized run candidates:

    ${groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(candidates))}
    """
}

workflow {
    // Work queue: a handful of real ENA run accessions chosen to exercise the messy
    // cases (edit this list to point at your own study). Drawn from two public RNA-seq
    // studies on purpose:
    //   SRR1039508  airway (SRP033351)  PAIRED, clean two-file record
    //   SRR1039513  airway (SRP033351)  PAIRED, ENA also lists an orphan .fastq.gz to drop
    //   SRR031708   pasilla (SRP001537) SINGLE  ) same biosample SAMN00006272
    //   SRR031712   pasilla (SRP001537) SINGLE  ) -> REDUCE merges as replicates
    //   SRR031714   pasilla (SRP001537) PAIRED
    //   SRR031726   pasilla (SRP001537) PAIRED  ) same biosample SAMN00006277
    //   SRR031727   pasilla (SRP001537) PAIRED  ) -> REDUCE merges as replicates
    def accessions = channel.of(
        'SRR1039508', 'SRR1039513',
        'SRR031708', 'SRR031712', 'SRR031714', 'SRR031726', 'SRR031727'
    )

    def meta  = ENA_FETCH(accessions)     // SHARD  (parallel process fan-out)
    def cands = normalize(meta)           // MAP    (one TaskRun per run, parallel up to executor.cpus)

    consolidate(cands.collect()).view { s ->   // REDUCE (gather -> one agent call)
        file("${projectDir}/samplesheet.csv").text = s.csv
        """\
        Wrote ${projectDir}/samplesheet.csv

        ${s.csv}
        --- notes ---
        ${s.notes}
        """.stripIndent()
    }
}
