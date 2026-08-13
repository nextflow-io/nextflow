nextflow.enable.types = true

// Fan-in parity: a canonical `process` and an `agent` consume the SAME collected
// channel (a `Bag<Finding>` produced by `collect()`) and BOTH fire exactly once over
// the whole bag. This shows the agent inherits canonical Nextflow cardinality: a
// value/singleton input (what `collect()` yields) => one invocation. The only
// difference is what the body does — deterministic Groovy vs an LLM call. See README.md.

record Finding {
    id:      String
    summary: String
}

// Canonical process reducer: one value input (the collected Bag) -> runs ONCE.
process reduce_with_process {
    input:
    findings: Bag<Finding>
    output:
    report: String
    exec:
        report = "combined ${findings.size()} findings [${findings.collect { it.id }.sort().join(', ')}]"
}

// Agent reducer: SAME input shape (Bag<Finding>) -> also runs ONCE (same isSingleton rule).
agent reduce_with_agent {
    model 'openai/gpt-5-mini'
    instruction 'You synthesise a set of findings into ONE short summary sentence.'
    input:
    findings: Bag<Finding>
    output:
    report: String
    prompt:
    """
    Combine these ${findings.size()} findings into one short sentence:

    ${findings.collect { "- ${it.id}: ${it.summary}" }.join('\n')}
    """
}

workflow {
    // A queue channel of 3 findings...
    def findings = channel.of(
        record(id: 'f1', summary: 'coverage looks good'),
        record(id: 'f2', summary: 'contamination is low'),
        record(id: 'f3', summary: 'adapter content is minimal')
    )

    // ...collected into ONE value item (a Bag of all three). A value/singleton channel
    // is broadcast to every consumer, so both reducers below read the same bag.
    def bag = findings.collect()

    reduce_with_process(bag).view { r -> "PROCESS reducer => ${r}" }
    reduce_with_agent(bag).view   { r -> "AGENT   reducer => ${r}" }
}
