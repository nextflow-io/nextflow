nextflow.enable.types = true

// Two no-tool agents chained over a channel: stage-1's output record type IS
// stage-2's input type, so structured output flows agent→agent. See README.md.

// The contract BETWEEN the two agents: stage-1 output == stage-2 input.
record Hypothesis {
    statement:  String
    rationale:  String
    confidence: Double      // the author agent's self-reported confidence, 0..1
}

// The final output of the pipeline.
record Review {
    verdict:  String        // e.g. "plausible", "doubtful", "needs revision"
    critique: String        // a concise peer-review note
}

// Stage 1: propose a hypothesis (String question -> Hypothesis record).
agent hypothesizer {
    model 'openai/gpt-5-mini'
    instruction 'You are a scientist. Given a question, propose ONE concrete, testable hypothesis.'

    input:
        question: String

    output:
        hypothesis: Hypothesis

    prompt:
    """
    Question: ${question}

    Propose a single, testable hypothesis. Include a brief rationale and your
    confidence (a number between 0 and 1).
    """
}

// Stage 2: peer-review the hypothesis (Hypothesis -> Review record).
agent critic {
    model 'openai/gpt-5-mini'
    instruction 'You are a skeptical peer reviewer. Judge a hypothesis for plausibility and testability.'

    input:
        hypothesis: Hypothesis

    output:
        review: Review

    prompt:
    """
    Hypothesis: ${hypothesis.statement}
    Rationale:  ${hypothesis.rationale}
    Author confidence: ${hypothesis.confidence}

    Give a one-word verdict and a concise critique (testability, confounders,
    whether the author's confidence is justified).
    """
}

workflow {
    // stage-1's output channel feeds straight into stage-2.
    def hypotheses = hypothesizer(channel.of(
        'Does caffeine improve short-term memory consolidation?',
        'Can a high-fiber diet reduce systemic inflammation?'
    ))

    critic(hypotheses).view { r ->
        """\
        verdict  : ${r.verdict}
        critique : ${r.critique}
        """.stripIndent()
    }
}
