nextflow.enable.types = true

/*
 * The deterministic half of the report: the verdict is decided by CODE, not by the
 * model. `nf:module_run` exposes this process to the agent as the `qc_verdict` tool
 * because the agent module includes it — no network. The agent task itself always runs
 * in the `pi` runner image, so the engine that image needs also runs this container.
 */
process qc_verdict {
    container 'ubuntu:24.04'
    input:
    n50_kb: Integer
    completeness_pct: Float

    output:
    verdict: String = stdout()

    script:
    """
    awk -v n50='${n50_kb}' -v comp='${completeness_pct}' 'BEGIN {
        if( n50 >= 40 && comp >= 95.0 )      v = "PASS"
        else if( n50 >= 20 && comp >= 90.0 ) v = "WARN"
        else                                 v = "FAIL"
        printf "%s", v
    }'
    """
}
