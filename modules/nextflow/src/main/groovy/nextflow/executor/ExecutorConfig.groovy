/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nextflow.executor

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.Duration
import nextflow.util.MemoryUnit

@ScopeName("executor")
@Description("""
    The `executor` scope controls various executor behaviors.
""")
@CompileStatic
class ExecutorConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        *Used only by the SLURM, LSF, PBS and PBS Pro executors.*

        The project or organization account that should be charged for running the pipeline jobs.
    """)
    final String account

    @ConfigOption
    @Description("""
        Determines how long to wait for a job array to fill to the specified array size before submitting a partial array (default: `60 min`).
    """)
    final Duration arrayTimeout

    @ConfigOption
    @Description("""
        *Used only by the local executor.*

        The maximum number of CPUs made available by the underlying system.
    """)
    final Integer cpus

    @ConfigOption
    @Description("""
        Determines how often to log the executor status (default: `5 min`).
    """)
    final Duration dumpInterval

    @ConfigOption
    @Description("""
        *Used only by grid executors.*

        Determines how long to wait for the `.exitcode` file to be created after the task has completed, before returning an error status (default: `270 sec`).
    """)
    final Duration exitReadTimeout

    @ConfigOption
    @Description('''
        *Used only by grid executors and Google Batch.*

        Determines the name of jobs submitted to the underlying cluster executor:
        ```groovy
        executor.jobName = { "$task.name - $task.hash" }
        ```
    ''')
    final Closure<String> jobName

    @ConfigOption
    @Description("""
        Determines the number of jobs that can be killed in a single command execution (default: `100`).
    """)
    final int killBatchSize

    @ConfigOption
    @Description("""
        *Used only by the local executor.*

        The maximum amount of memory made available by the underlying system.
    """)
    final MemoryUnit memory

    @ConfigOption
    @Description("""
        The name of the executor to be used (default: `local`).
    """)
    final String name

    @ConfigOption
    @Description("""
        *Used only by the SLURM executor.*

        When `true`, memory allocations for SLURM jobs are specified as `--mem-per-cpu <task.memory / task.cpus>` instead of `--mem <task.memory>`.
    """)
    final boolean perCpuMemAllocation

    @ConfigOption
    @Description("""
        *Used only by the LSF executor.*

        Enables the *per-job* memory limit mode for LSF jobs.
    """)
    final boolean perJobMemLimit

    @ConfigOption
    @Description("""
        *Used only by the SLURM executor.*

        When `true`, job status queries use `squeue --only-job-state` without partition (`-p`) or user (`-u`) filters. This is required for SLURM installations where `squeue` only accepts the `--only-job-state` option (default: `false`).
    """)
    final boolean onlyJobState

    @ConfigOption
    @Description("""
        *Used only by the LSF executor.*

        Enables the *per-task* memory reserve mode for LSF jobs.
    """)
    final boolean perTaskReserve

    @ConfigOption
    @Description("""
        Determines how often to check for process termination. Default varies for each executor.
    """)
    final Duration pollInterval

    @ConfigOption
    @Description("""
        Determines how job status is retrieved. When `false` only the queue associated with the job execution is queried. When `true` the job status is queried globally i.e. irrespective of the submission queue (default: `false`).
    """)
    final boolean queueGlobalStatus

    @ConfigOption
    @Description("""
        The number of tasks the executor will handle in a parallel manner. A queue size of zero corresponds to no limit. Default varies for each executor.
    """)
    final Integer queueSize

    @ConfigOption
    @Description("""
        *Used only by grid executors.*

        Determines how often to fetch the queue status from the scheduler (default: `1 min`).
    """)
    final Duration queueStatInterval
    
    @Description("""
        The `executor.retry` scope controls the behavior of retrying failed job submissions.

        [Read more](https://nextflow.io/docs/latest/reference/config.html#executor)
    """)
    final ExecutorRetryConfig retry

    @ConfigOption
    @Description("""
        Determines the max rate of job submission per time unit, for example `'10sec'` (10 jobs per second) or `'50/2min'` (50 jobs every 2 minutes) (default: unlimited).
    """)
    final String submitRateLimit

    private Map opts

    /* required by extension point -- do not remove */
    ExecutorConfig() {}

    ExecutorConfig(Map opts) {
        account = opts.account
        arrayTimeout = opts.arrayTimeout as Duration ?: Duration.of('60min')
        cpus = opts.cpus as Integer
        dumpInterval = opts.dumpInterval as Duration ?: Duration.of('5min')
        exitReadTimeout = opts.exitReadTimeout as Duration ?: Duration.of('270sec')
        jobName = opts.jobName as Closure
        killBatchSize = opts.killBatchSize != null ? opts.killBatchSize as int : 100
        memory = opts.memory as MemoryUnit
        name = opts.name
        perCpuMemAllocation = opts.perCpuMemAllocation as boolean
        perJobMemLimit = opts.perJobMemLimit as boolean
        onlyJobState = opts.onlyJobState as boolean
        perTaskReserve = opts.perTaskReserve as boolean
        pollInterval = opts.pollInterval as Duration
        queueGlobalStatus = opts.queueGlobalStatus as boolean
        queueSize = opts.queueSize as Integer
        queueStatInterval = opts.queueStatInterval as Duration ?: Duration.of('1min')
        retry = opts.retry ? new ExecutorRetryConfig(opts.retry as Map) : new ExecutorRetryConfig(Map.of())
        submitRateLimit = opts.submitRateLimit

        // preserve executor-specific opts
        this.opts = opts
    }

    Duration getArrayTimeout(String execName) {
        getExecConfigProp(execName, 'arrayTimeout', null) as Duration
    }

    Duration getExitReadTimeout(String execName) {
        getExecConfigProp(execName, 'exitReadTimeout', null) as Duration
    }

    Duration getMonitorDumpInterval(String execName) {
        getExecConfigProp(execName, 'dumpInterval', null) as Duration
    }

    Duration getPollInterval(String execName, Duration defValue = null) {
        getExecConfigProp(execName, 'pollInterval', defValue ?: Duration.of('1sec')) as Duration
    }

    int getQueueSize(String execName, int defValue) {
        getExecConfigProp(execName, 'queueSize', defValue) as int
    }

    Duration getQueueStatInterval(String execName) {
        getExecConfigProp(execName, 'queueStatInterval', null) as Duration
    }

    @Memoized
    Object getExecConfigProp(String execName, String name, Object defValue, Map env = null) {
        // -- check `executor.$<execName>.<name>`
        final execProp = execProp(execName, name)
        if( execProp != null )
            return execProp

        // -- check `executor.<name>`
        final prop = this.hasProperty(name) ? this.getProperty(name) : null
        if( prop != null )
            return prop

        // -- check environment variable
        final key = "NXF_EXECUTOR_${name.toUpperCase().replaceAll(/\./,'_')}".toString()
        if( env == null )
            env = System.getenv()
        return env.containsKey(key) ? env.get(key) : defValue
    }

    private Object execProp(String execName, String name) {
        if( !execName )
            return null
        final result = opts['$' + execName]
        return result instanceof Map ? result[name] : null
    }

}
