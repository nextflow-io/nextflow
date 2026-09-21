nextflow.enable.types = true

// Fully-agentic map-reduce: every phase is an `agent`, composed with plain channels.
//
//   PLAN   (agent)  : planner breaks a brief into independent research shards (one LLM call).
//   SHARD  (channel): plan.plan.flatMap { it.shards } — deterministic work queue.
//   MAP    (agent)  : mapper answers each shard — one TaskRun per shard, run in PARALLEL
//                     (bounded by maxForks / executor.local.cpus).
//   REDUCE (agent)  : reducer synthesises one report from ALL findings via collect() (fan-in).
//
// This is the shape from the process-parity design (spec §13): the agent primitive now runs as
// a real task, so the map fans out concurrently, each call is cached/resumable, and every box
// shows up in the progress table and lineage. Contrast `samplesheet-builder`, where the SHARD is
// a deterministic process; here the plan itself is agentic. See README.md.

record Shard   { id: String;      question: String }
record Plan    { title: String;   shards: List<Shard> }
record Finding { shardId: String; summary: String }
record Report  { title: String;   body: String }

// PLAN — one value input -> runs once -> planner.out.plan (single record, bare schema).
agent planner {
    // NOTE: a floating alias. The Pi model registry currently ships no dated `gpt-5*`
    // snapshot, and every dated `gpt-4o*` id it does know returns HTTP 500 from OpenAI's
    // Responses API (the only API the Pi runner uses). Pin a dated snapshot here once one
    // is available: the cache key follows the alias, so a server-side model upgrade makes
    // `-resume` replay a generation the current model would not produce.
    model 'openai/gpt-5-mini'
    instruction '''\
        You decompose a research brief into a small set of INDEPENDENT sub-questions
        ("shards") that can be investigated in parallel without depending on each other.
        Give each shard a short stable id and one focused question. Prefer 3-5 shards.
        '''.stripIndent()

    input:
    brief: String
    output:
    plan: Plan

    prompt:
    """
    Break this brief into independent research shards:

    ${brief}
    """
}

// MAP — one shard per invocation -> one TaskRun each (parallel).
agent mapper {
    model 'openai/gpt-5-mini'
    instruction 'You are a focused researcher. Answer the single question concisely and factually, and echo the shard id you were given so results can be correlated.'

    input:
    shard: Shard
    output:
    finding: Finding

    prompt:
    """
    Shard id: ${shard.id}
    Answer this shard question:

    ${shard.question}
    """
}

// REDUCE — the whole bag of findings as ONE value input -> fires once (fan-in).
agent reducer {
    model 'openai/gpt-5-mini'
    instruction 'You synthesise many independent findings into one coherent, non-redundant report.'

    input:
    findings: Bag<Finding>     // `collect()` yields a Bag (order-independent) — fan-in of all findings
    output:
    report: Report

    prompt:
    """
    Synthesise ONE report from these findings:

    ${findings.collect { "- (${it.shardId}) ${it.summary}" }.join('\n')}
    """
}

workflow {
    def brief = params.brief   // default in nextflow.config; override with --brief '...'

    // A single-output agent auto-unwraps to its output channel under the typed DSL,
    // so `planner(...)` IS the channel of Plan records (no `.plan` accessor needed).
    def plan     = planner(channel.of(brief))       // PLAN   (agentic -> channel<Plan>)
    def shards   = plan.flatMap { it.shards }        // SHARD  (deterministic work queue -> channel<Shard>)
    def findings = mapper(shards)                    // MAP    (one TaskRun per shard, parallel -> channel<Finding>)

    reducer(findings.collect())                      // REDUCE (collect() -> one value -> single TaskRun)
        .view { r ->
            """\
            ===== ${r.title} =====
            ${r.body}
            """.stripIndent()
        }
}
