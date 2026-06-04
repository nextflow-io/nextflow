nextflow.enable.types = true

/*
 * Input record. The agent receives one `Query` per channel item. The optional
 * `context` field (note the `?`) is honored on the input side.
 */
record Query {
    question: String
    context: String?
}

/*
 * Output record. Its structure is reflected into a JSON schema that constrains
 * the model's response (OpenAI structured output). The returned JSON is parsed
 * and bound to an `Analysis` record emitted on the channel. The field types
 * below exercise the full v1 schema deriver: String, Double, Boolean and a
 * List of strings.
 */
record Analysis {
    summary: String
    confidence: Double
    actionable: Boolean
    key_points: List<String>
}

agent analyst {
    model 'openai/gpt-5-mini'
    instruction 'You are a precise scientific analyst. Be concise and honest about uncertainty.'

    input:
        query: Query

    output:
        result: Analysis

    prompt:
    """
    Analyze the following question and return a structured analysis.

    Question: ${query.question}
    """
}

workflow {
    // Typed DSL: invoke with the call form and chain operators with `.view { }`
    // (the bare `| view` pipe is not used under `nextflow.enable.types = true`).
    analyst(channel.of(
        record(question: 'Is FASTQ a binary or a text format?', context: 'bioinformatics file formats')
    ))
    .view { r ->
        """\
        summary    : ${r.summary}
        confidence : ${r.confidence}
        actionable : ${r.actionable}
        key_points : ${r.key_points}
        """.stripIndent()
    }
}
