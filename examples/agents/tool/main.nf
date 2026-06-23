nextflow.enable.types = true

/*
 * A normal Nextflow process. The agent exposes it to the LLM as a tool via
 * the `module_run` capability: the LLM's tool-call args are marshalled into
 * this process's input channel, the process runs as a real dataflow node
 * (executor / work dir / cache), and its output is serialized back to the LLM
 * as JSON.
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
 * The agent. `tools 'module_run'` enables the generic `module_run(module, args)`
 * tool. Because `uppercase` is an in-scope process (defined above in the same
 * script), it is auto-discovered — no `include` statement needed. The LLM
 * receives a `module_run` tool with a `module` enum listing `uppercase`.
 *
 * The LLM decides when to call module_run; the harness runs the process as a
 * dataflow node and feeds the result back so the LLM can produce its final
 * answer.
 */
agent shouty {
    model 'openai/gpt-5-mini'
    instruction 'To uppercase text you MUST call module_run with {"module":"uppercase","args":{"text":"<input>"}}, then reply with only the result.'
    tools 'module_run'

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
