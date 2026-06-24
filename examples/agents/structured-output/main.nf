nextflow.enable.types = true

// Input record — `context` is optional (note `?`).
record Query {
    question: String
    context: String?
}

// Output record — its fields become the model's JSON-schema contract (structured output).
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
    // Under typed DSL, use the call form; `| view` pipe is not available.
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
