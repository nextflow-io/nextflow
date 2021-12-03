/*
 * Copyright 2020-2021, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
import java.nio.file.Path
import java.util.regex.Pattern

import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
/**
 * Processor for SLURM resource manager (DRAFT)
 *
 * See http://computing.llnl.gov/linux/slurm/
 *
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class SlurmExecutor extends AbstractGridExecutor {

    static private Pattern SUBMIT_REGEX = ~/Submitted batch job (\d+)\s*$/
    static private Pattern SUBMIT_REGEX_WITH_CLUSTER = ~/Submitted batch job (\d+) on cluster (\S+)\s*$/

    private boolean hasSignalOpt(Map config) {
        def opts = config.clusterOptions?.toString()
        return opts ? opts.contains('--signal ') || opts.contains('--signal=') : false
    }


    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        result << '-D' << quote(task.workDir)
        result << '-J' << getJobNameFor(task)
        result << '-o' << quote(task.workDir.resolve(TaskRun.CMD_LOG))     // -o OUTFILE and no -e option => stdout and stderr merged to stdout/OUTFILE
        result << '--no-requeue' << '' // note: directive need to be returned as pairs

        if( !hasSignalOpt(task.config) ) {
            // see https://github.com/nextflow-io/nextflow/issues/2163
            // and https://slurm.schedmd.com/sbatch.html#OPT_signal
            result << '--signal' << 'B:USR2@30'
        }

        if( task.config.cpus > 1 ) {
            result << '-c' << task.config.cpus.toString()
        }

        if( task.config.time ) {
            result << '-t' << task.config.getTime().format('HH:mm:ss')
        }

        if( task.config.getMemory() ) {
            //NOTE: Enforcement of memory limits currently relies upon the task/cgroup plugin or
            // enabling of accounting, which samples memory use on a periodic basis (data need not
            // be stored, just collected). In both cases memory use is based upon the job's
            // Resident Set Size (RSS). A task may exceed the memory limit until the next periodic
            // accounting sample. -- https://slurm.schedmd.com/sbatch.html
            result << '--mem' << task.config.getMemory().toMega().toString() + 'M'
        }

        // the requested partition (a.k.a queue) name
        if( task.config.queue ) {
            result << '-p' << (task.config.queue.toString())
        }

        // the requested cluster name
        if( task.config.cluster ) {
            result << '-M' << (task.config.cluster.toString())
        }

        // -- at the end append the command script wrapped file name
        if( task.config.clusterOptions ) {
            result << task.config.clusterOptions.toString() << ''
        }

        return result
    }

    String getHeaderToken() { '#SBATCH' }

    /**
     * The command line to submit this job
     *
     * @param task The {@link TaskRun} instance to submit for execution to the cluster
     * @param scriptFile The file containing the job launcher script
     * @return A list representing the submit command line
     */
    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {

        ['sbatch', scriptFile.getName()]

    }

    /**
     * Parse the string returned by the {@code sbatch} command and extract the job ID string
     *
     * @param text The string returned when submitting the job
     * @return The actual job ID string
     */
    @Override
    def parseJobId(String text) {

        for( String line : text.readLines() ) {
            def m = SUBMIT_REGEX.matcher(line)
            def m2 = SUBMIT_REGEX_WITH_CLUSTER.matcher(line)
            if( m.find() ) {
                return m.group(1).toString()
            }
            if( m2.find() ) {
                return new JobIdOnCluster(jobId: m2.group(1).toString(), cluster: m2.group(2).toString())
            }
        }

        // customised `sbatch` command can return only the jobid
        def id = text.trim()
        if( id.isLong() )
            return id
        def idc = id.split(/;/) // This is returned when passing --parsable to sbatch
        if(idc.size() == 2 && idc[0].isLong())
            return new JobIdOnCluster(jobId: idc[0], cluster: idc[1])
        throw new IllegalStateException("Invalid SLURM submit response:\n$text\n\n")
    }

    def killCommandForCluster(Optional<String> cluster, List<Object> jobId) {
        final result = getKillCommand()
        if(cluster.isPresent()) {
            result.add('-M')
            result.add(cluster.get())
        }
        result.addAll(jobId.collect { it.toString()})
        return result
    }

    def execAndGetReturnValue(List<String> cmd) {
        def proc = new ProcessBuilder(cmd).redirectErrorStream(true).start()
        proc.waitForOrKill(10_000)
        return proc
    }

    @Override
    void killTask(Object jobId) {
        List<List<String>> execList = getExecListForKill(jobId)
        execList.collect {cmd ->
            [cmd: cmd, proc: execAndGetReturnValue(cmd)]
        }.findAll {
            it.get('proc').exitValue() != 0
        }.each {it ->
            reportFailure(cmd: it.get('cmd'), proc: it.get('proc'))
        }
    }

    protected List<List<String>> getExecListForKill(jobId) {
        def jobIds = jobId instanceof Collection ? jobId : [jobId]
        def richJobIds = jobIds.findAll { it instanceof JobIdOnCluster }
        def byCluster = richJobIds.inject([:]) {
            m, it -> m << [(it.cluster):  m.getOrCreate(it.cluster, []) << it.jobId]
        }
        def execList = byCluster.collect {
            c, jids -> killCommandForCluster(Optional.of(c), jids)
        }
        def simpleJobIds = jobIds.findAll { !(it instanceof JobIdOnCluster) }
        if(simpleJobIds) execList.add(killCommandForCluster(Optional.empty(), simpleJobIds))
        execList
    }

    @Override
    protected List<String> getKillCommand() { ['scancel'] }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        final clusterSuffix = queue && queue instanceof QueueOnCluster ? ";${queue.cluster}".toString() : ''

        final result = ['squeue','--noheader','-o',"%i${clusterSuffix} %t".toString(), '-t', 'all']

        if( queue )
            if (queue instanceof QueueOnCluster) {
                result << '-p' << queue.queueName
                result << '-M' << queue.cluster
            }
            else
                result << '-p' << queue.toString()

        final user = System.getProperty('user.name')
        if( user )
            result << '-u' << user
        else
            log.debug "Cannot retrieve current user"

        return result
    }

    /*
     *  Maps SLURM job status to nextflow status
     *  see http://slurm.schedmd.com/squeue.html#SECTION_JOB-STATE-CODES
     */
    static private Map STATUS_MAP = [
            'PD': QueueStatus.PENDING,  // (pending)
            'R': QueueStatus.RUNNING,   // (running)
            'CA': QueueStatus.ERROR,    // (cancelled)
            'CF': QueueStatus.PENDING,  // (configuring)
            'CG': QueueStatus.RUNNING,  // (completing)
            'CD': QueueStatus.DONE,     // (completed)
            'F': QueueStatus.ERROR,     // (failed),
            'TO': QueueStatus.ERROR,    // (timeout),
            'NF': QueueStatus.ERROR,    // (node failure)
            'S': QueueStatus.HOLD,      // (job suspended)
            'ST': QueueStatus.HOLD,     // (stopped)
            'PR': QueueStatus.ERROR,    // (Job terminated due to preemption)
            'BF': QueueStatus.ERROR,    // (boot fail, Job terminated due to launch failure)
    ]

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {

        def result = [:]

        text.eachLine { String line ->
            def cols = line.split(/\s+/)
            if( cols.size() == 2 ) {
                def jobId = cols[0].split(/;/)
                if (jobId.size() == 2) {
                    def jidc = new JobIdOnCluster(jobId: jobId[0], cluster: jobId[1])
                    result.put(jidc.toString(), STATUS_MAP.get(cols[1]) )
                }
                else {
                    result.put( cols[0], STATUS_MAP.get(cols[1]) )
                }
            }
            else {
                log.debug "[SLURM] invalid status line: `$line`"
            }
        }

        return result
    }

    void reportFailure(List<String> cmd, Process proc) {
        def ret = proc.exitValue()
        def m = """\
                Unable to kill pending jobs
                - cmd executed: ${cmd.join(' ')}
                - exit status : $ret
                - output      :
                """
                .stripIndent()
        m += proc.text.indent('  ')
        log.debug(m)
    }
}

/**
 * Specialises queue with cluster
 */
@Slf4j
class QueueOnCluster {
    String queueName
    String cluster
    @Override
    String toString() {"${queueName}:${cluster}"}
}

/**
 * Specialises JobId with cluster
 */
@Slf4j
class JobIdOnCluster {
    Object jobId
    String cluster

    @Override
    String toString() {"${jobId};${cluster}"}

    @Override
    boolean equals(Object other) {
        if(other instanceof JobIdOnCluster) {
            return other.jobId == jobId && other.cluster == cluster
        }
        return false
    }

    @Override
    int hashCode() {
        return toString().hashCode()
    }
}