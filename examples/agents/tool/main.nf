nextflow.enable.types = true

/*
 * A normal Nextflow process. The agent exposes it to the LLM as its OWN tool
 * (named `uppercase`) via the `module_run` capability: the tool's parameters
 * schema IS this process's input schema ({text}, required), so the LLM's call
 * is validated against it. The args are marshalled into this process's input
 * channel, the process runs as a real dataflow node (executor / work dir /
 * cache), and its output is serialized back to the LLM as JSON.
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
 * The agent. `tools 'module_run'` exposes every in-scope process as its OWN
 * tool, named after the module. Because `uppercase` is an in-scope process
 * (defined above in the same script), it is auto-discovered — no `include`
 * statement needed. The LLM receives a tool named `uppercase` whose parameters
 * schema is {text} (required).
 *
 * The LLM decides when to call the `uppercase` tool; the harness runs the
 * process as a dataflow node and feeds the result back so the LLM can produce
 * its final answer.
 */
agent shouty {
    model 'openai/gpt-5-mini'
    instruction 'To uppercase text you MUST call the `uppercase` tool with {"text":"<input>"}, then reply with only the result.'
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
