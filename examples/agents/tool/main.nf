nextflow.enable.types = true

// A local process; `module_run` exposes it to the LLM as the tool `uppercase`.
process uppercase {
    input:
        text: String
    output:
        result: String
    exec:
        result = text.toUpperCase()
}

// `tools 'module_run'` auto-discovers in-scope processes and exposes each as a tool.
agent shouty {
    model 'openai/gpt-5-mini'
    instruction 'To uppercase text, call the `uppercase` tool, then reply with only the result.'
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
