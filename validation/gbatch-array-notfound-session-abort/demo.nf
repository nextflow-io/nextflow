/*
 * Reproduces the bug fixed by this PR: a Google Batch job array aborts the WHOLE
 * Nextflow session when one of its tasks exceeds `maxSubmitAwait`, instead of just
 * retrying that one task.
 *
 * All 3 tasks below share ONE Google Batch job, because of `array 3`. The first
 * submit attempt is given only 5 seconds to start -- far less than a VM needs to
 * provision -- so it deterministically times out (no real capacity shortage needed).
 * Nextflow then deletes the shared job and re-submits the timed-out task (attempt 2,
 * a generous 30 minute limit).
 *
 * Before this fix: the other tasks in the array poll the job Nextflow just deleted,
 * get NOT_FOUND, and that exception escapes uncaught -- aborting the whole run and
 * killing every other in-flight task, not just this array. Look for "Session
 * aborted" in the log.
 *
 * After this fix: the run finishes normally (exit code 0, "Execution complete --
 * Goodbye" in the log, no "Session aborted"). You'll still see a FAILED trace record
 * per task for its timed-out first attempt -- that's the normal, expected retry
 * bookkeeping, not the bug; what matters is that the run itself doesn't abort and all
 * 3 "task N ran to completion" lines print.
 */
process demo {
    array 3
    container 'ubuntu:24.04'
    errorStrategy 'retry'
    maxRetries 3
    maxSubmitAwait { task.submitAttempt == 1 ? 5.s : 30.m }

    input:
    val i

    output:
    stdout

    script:
    """
    echo "task ${i} ran to completion"
    """
}

workflow {
    channel.of(1..3) | demo | view
}
