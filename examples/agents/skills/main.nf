nextflow.enable.types = true

/*
 * Agent skills example (Milestone 4).
 *
 * The `skills` directive gives the agent access to one or more Anthropic-style
 * skills (SKILL.md folders). This agent declares the local `sequence-report`
 * skill found under `skills/sequence-report/`. When the prompt calls for an
 * assembly/QC summary the model activates the skill (langchain4j Tool Mode:
 * `activate_skill`) and follows its instructions — producing the distinctive
 * `[SEQ-REPORT v1]` formatted report it would not otherwise emit.
 *
 * A skill entry may also be a remote GitHub reference, e.g.
 *     skills 'github.com/<org>/<repo>@<sha>'
 * which is cloned and cached into this same `skills/` directory.
 */
agent reporter {
    model 'openai/gpt-5-mini'
    instruction 'You summarize sequencing and genome-assembly results for bioinformaticians.'
    skills 'sequence-report'

    input:
        request: String
    output:
        answer: String

    prompt:
    """
    ${request}
    """
}

workflow {
    reporter(channel.of('Summarize this assembly: N50 = 45 kb, total length = 5.1 Mb, GC = 50.8%. Is it good enough to proceed?'))
        .view { a -> "ANSWER=\n${a}" }
}
