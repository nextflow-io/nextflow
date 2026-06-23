nextflow.enable.types = true

/*
 * ============================================================================
 *  TWO AGENTS, NO TOOLS, CONNECTED BY A CHANNEL
 * ============================================================================
 *
 *  Agents compose with the rest of a pipeline through the normal channel /
 *  workflow model -- exactly like processes. This example chains two NO-TOOL
 *  agents: the first proposes a hypothesis, the second peer-reviews it. The
 *  output of `hypothesizer` is wired straight into `critic` as its input.
 *
 *  Neither agent declares `tools`, so each runs a single-shot LLM call (no
 *  tool-call loop). Both use STRUCTURED output: the output record type is
 *  reflected into a JSON schema the model is constrained to, and the parsed
 *  JSON is bound to a typed record on the output channel.
 *
 *  The key idea: the FIRST agent's output record type (`Hypothesis`) is also
 *  the SECOND agent's input record type. That shared contract is what lets one
 *  agent feed the other over a channel -- no tools, no glue code.
 */

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

/*
 * Stage 1 -- propose a hypothesis. Input is a plain question String; output is
 * a structured `Hypothesis` record.
 */
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

/*
 * Stage 2 -- peer-review the hypothesis. Input is the `Hypothesis` record
 * emitted by stage 1; output is a structured `Review` record.
 */
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
    // Typed DSL: use the call form (the bare `| pipe` is not used under
    // `nextflow.enable.types = true`). The ChannelOut returned by the first
    // agent is passed straight into the second -- the `Hypothesis` records
    // flow from `hypothesizer` to `critic` over the channel.
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
