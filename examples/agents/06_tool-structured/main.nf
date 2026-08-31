nextflow.enable.types = true

// A local process; `nf:module_run` exposes it to the LLM as the tool `uppercase`.
process uppercase {
    input:
    text: String
    output:
    result: String
    exec:
        result = text.toUpperCase()
}

// A record output type: combining `tools` with a structured (record) output is
// supported (M5). The tool loop runs schema-free; a final structuring turn then
// converts the loop's free-text answer into schema-valid JSON that binds to the
// record. NOTE: the structuring turn is +1 LLM call and is potentially lossy
// (it re-encodes free text to JSON), so keep the output schema close to the
// tool's own result.
record Shout {
    result: String
}

agent shouty {
    model 'openai/gpt-5-mini'
    instruction 'To uppercase text, call the `uppercase` tool, then reply with the result.'
    tools 'nf:module_run'

    input:
    request: String
    output:
    answer: Shout

    prompt:
    """
    ${request}
    """
}

workflow {
    shouty(channel.of('uppercase the word hello'))
        .view { a -> "ANSWER=${a}" }
}
