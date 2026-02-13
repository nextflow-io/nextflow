/*
 * Test: Azure Batch onAllTasksComplete race condition
 *
 * Each Nextflow *process* maps to one Azure Batch *job*.
 * Each *task instance* of that process is an Azure Batch *task* within the job.
 * Pool task slots = number of CPUs on the node (not configurable independently).
 *
 * Scenario 1 (fast_tasks): Many trivially fast task instances in one job.
 *   Race window: task 1 completes before Nextflow submits task 2.
 *   Azure sees "1/1 tasks completed" → terminates job → task 2 gets 409.
 *   Expected: recreateJobForTask handles this via job recreation.
 *
 * Scenario 2 (slow_tasks): More task instances than available pool slots.
 *   Tasks that can't run immediately sit in 'active' state (queued).
 *   Tests Ghislain's claim: does Azure terminate the job while tasks
 *   are still in active/running state (not completed)?
 *   According to docs, onAllTasksComplete fires only when ALL tasks
 *   reach 'completed' state, so this should NOT cause premature termination.
 *
 * Usage:
 *   ./launch.sh run tests/azure-terminate-race.nf \
 *     -profile azure \
 *     -w az://yourcontainer/work \
 *     --num_fast 50 \
 *     --num_slow 20 \
 *     --slow_seconds 60
 *
 * To force queuing in Scenario 2, ensure num_slow exceeds
 * total task slots in the pool (nodes × CPUs per node).
 * e.g. pool with 1× Standard_D2s_v3 (2 CPUs) = 2 slots → use --num_slow 20
 *
 * Enable debug logging to see 409s and job recreation:
 *   NXF_DEBUG=2 ./launch.sh run tests/azure-terminate-race.nf ...
 */

params.num_fast = 50
params.num_slow = 20
params.slow_seconds = 60

/*
 * Scenario 1: Fast tasks — triggers the "all tasks completed" race
 *
 * Each instance is one Azure Batch task in a single job.
 * echo completes almost instantly, so task N likely finishes
 * before Nextflow submits task N+1. With TERMINATE_JOB set
 * eagerly, the job terminates and the next submission gets 409.
 */
process fast_tasks {
    cpus 1
    tag "${id}"

    input:
    val id

    output:
    stdout

    script:
    """
    echo "fast task ${id} done"
    """
}

/*
 * Scenario 2: Slow tasks — tests active-state termination claim
 *
 * Submits more tasks than pool slots. Excess tasks queue in
 * 'active' state. If Azure terminates the job while active
 * tasks exist, the pipeline will fail with missing outputs.
 *
 * cpus 1 ensures each task occupies exactly 1 slot, so a
 * 2-CPU node runs 2 tasks in parallel and queues the rest.
 */
process slow_tasks {
    cpus 1
    tag "${id}"

    input:
    val id

    output:
    stdout

    script:
    """
    echo "slow task ${id} started at \$(date)"
    sleep ${params.slow_seconds}
    echo "slow task ${id} finished at \$(date)"
    """
}

workflow {
    // Scenario 1: rapid-fire trivial tasks (same process = same job)
    fast_ch = Channel.of(1..params.num_fast)
    fast_results = fast_tasks(fast_ch)
    fast_results.count().view { n -> "Scenario 1 PASSED: ${n}/${params.num_fast} fast tasks completed" }

    // Scenario 2: slow tasks that exceed pool capacity (same process = same job)
    slow_ch = Channel.of(1..params.num_slow)
    slow_results = slow_tasks(slow_ch)
    slow_results.count().view { n -> "Scenario 2 PASSED: ${n}/${params.num_slow} slow tasks completed" }
}
