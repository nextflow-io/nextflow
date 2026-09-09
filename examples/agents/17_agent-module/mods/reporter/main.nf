/*
 * A self-contained agent module.
 *
 * Everything the agent needs is inside THIS directory:
 *   skills/qa-report/SKILL.md   report format
 *   skills/style/SKILL.md       house style
 *   tools/qc_verdict.nf         the deterministic PASS/WARN/FAIL rule
 *
 * Both are resolved relative to the module directory — the directory of the file
 * that DECLARES the agent — not the directory of whoever includes it, and that
 * holds under an alias too. So this module keeps behaving the same wherever it
 * is dropped, and a `skills/qa-report/` sitting next to the consumer's script
 * cannot shadow the module's own.
 *
 * `tools 'nf:module_run'` sees only what THIS file defines or includes, which is why
 * the module includes its own tool. A process defined by the including script is
 * deliberately invisible.
 *
 * No `nextflow.enable.types = true` here: an `agent` block does not need it. The
 * tool module does, because it declares a typed process — the flag is per file.
 */

include { qc_verdict } from './tools/qc_verdict.nf'

agent reporter {
    model 'openai/gpt-5-mini'
    instruction """
        You report on genome assembly QC for bioinformaticians.
        The PASS/WARN/FAIL verdict is NOT yours to guess: always obtain it by
        calling the qc_verdict tool with the metrics you were given.
        """

    skills 'qa-report', 'style'
    tools 'nf:module_run'

    input:
    request: String

    output:
    answer: String

    prompt:
    """
    ${request}
    """
}
