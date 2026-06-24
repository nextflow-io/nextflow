nextflow.enable.types = true

// Convergence loop: the agent re-runs the SAME tool (score_threshold) many times,
// varying `threshold`, reading the F1 each call, converging on the optimum. The knob
// is a declared input; the dataset path is fixed (params.scores) so each tool call is
// just {"threshold": <x>}. Fully offline. See README.md.

// Tool `score_threshold`: precision/recall/F1 of `score >= threshold` on the dataset.
process score_threshold {
    input:
        threshold: BigDecimal
    output:
        result: String
    exec:
        def tp = 0
        def fp = 0
        def fn = 0
        def t = threshold as Double
        new File(params.scores as String).eachLine { line, n ->
            if( n == 1 ) return                       // skip header
            def parts = line.trim().split(/\s+/)
            if( parts.size() < 2 ) return
            def score = parts[0] as Double
            def label = parts[1] as Integer
            def pred = score >= t ? 1 : 0
            if( pred == 1 && label == 1 ) { tp = tp + 1 }
            if( pred == 1 && label == 0 ) { fp = fp + 1 }
            if( pred == 0 && label == 1 ) { fn = fn + 1 }
        }
        def prec = (tp + fp) > 0 ? (tp / (tp + fp)) : 0.0
        def rec  = (tp + fn) > 0 ? (tp / (tp + fn)) : 0.0
        def f1   = (prec + rec) > 0 ? (2 * prec * rec / (prec + rec)) : 0.0
        result = String.format('threshold=%.3f tp=%d fp=%d fn=%d precision=%.4f recall=%.4f f1=%.4f', t, tp, fp, fn, prec as Double, rec as Double, f1 as Double)
}

agent tuner {
    model 'openai/gpt-5-mini'

    // role/constraints only — no tool-field structure, no fixed step list.
    instruction '''\
        You optimise a decision threshold. The score_threshold tool is the ONLY
        way to learn the F1 for a given threshold — never guess F1, always call
        the tool. Search efficiently: scan the range coarsely, then refine around
        the best value until further changes do not improve F1.
        '''.stripIndent()

    // the OBJECTIVE that drives the convergence loop.
    goal '''\
        Find the score threshold in [0, 1] that MAXIMISES F1. Call score_threshold
        with different thresholds, read the F1 it returns, and converge on the
        best one. Report the optimal threshold and its F1 (with precision and
        recall).
        '''.stripIndent()

    tools 'module_run'
    maxIterations 30

    input:
        request: String
    output:
        report: String

    prompt:
    """
    ${request}
    """
}

workflow {
    // The dataset (params.scores, default below) is committed alongside this
    // example — small, self-contained, no container or network needed.
    tuner(channel.of('Find the score threshold that maximises F1 on the labelled dataset.'))
        .view { r -> "OPTIMAL=${r}" }
}
