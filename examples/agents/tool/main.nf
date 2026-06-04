nextflow.enable.types = true

/*
 * A normal Nextflow process. The agent exposes it to the LLM as a tool:
 * the LLM's tool-call args are marshalled into this process's input channel,
 * the process runs as a real dataflow node (executor / work dir / cache),
 * and its output is serialized back to the LLM as JSON.
 */
process uppercase {
    input:
        text: String
    output:
        result: String
    exec:
        result = text.toUpperCase()
}

/*
 * The agent. `tools 'uppercase'` references the in-scope process above.
 * The LLM decides when to call it; the harness runs it and feeds the
 * result back so the LLM can produce its final answer.
 */
agent shouty {
    model 'openai/gpt-5-mini'
    instruction 'To uppercase text you MUST call the uppercase tool, then reply with only its result.'
    tools 'uppercase'

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
    shouty(channel.of('uppercase the word hello'))
        .view { a -> "ANSWER=${a}" }
}
