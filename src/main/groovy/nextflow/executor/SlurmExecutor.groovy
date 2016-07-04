/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.executor
import java.nio.file.Path

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

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        result << '-D' << task.workDir.toString()
        result << '-J' << getJobNameFor(task)
        result << '-o' << task.workDir.resolve(TaskRun.CMD_LOG).toString()           // -o OUTFILE and no -e option => stdout and stderr merged to stdout/OUTFILE

        if( task.config.cpus > 1 ) {
            result << '-c' << task.config.cpus.toString()
        }

        if( task.config.time ) {
            result << '-t' << task.config.getTime().format('HH:mm:ss')
        }

        if( task.config.getMemory() ) {
            result << '--mem' << task.config.getMemory().toMega().toString()
        }

        // the requested partition (a.k.a queue) name
        if( task.config.queue ) {
            result << '-p' << (task.config.queue.toString())
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
        def pattern = ~ /Submitted batch job (\d+)/
        for( String line : text.readLines() ) {
            def m = pattern.matcher(line)
            if( m.matches() ) {
                return m.group(1).toString()
            }
        }

        throw new IllegalStateException("Invalid SLURM submit response:\n$text\n\n")
    }

    protected String getKillCommand() { 'scancel' }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        if( queue )
            log.debug "SLURM executor does not support queue parameter on queue status"

        return ['squeue','-h','-o','%i %t', '-t', 'all']
    }

    /*
     *  Maps SLURM job status to nextflow status
     *  see https://computing.llnl.gov/linux/slurm/squeue.html#SECTION_JOB%20STATE%20CODES
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
            'PR': QueueStatus.ERROR,    // (Job terminated due to preemption)
    ]

    @Override
    protected Map<?, QueueStatus> parseQueueStatus(String text) {

        def result = [:]

        text.eachLine { String line ->
            def cols = line.split(/\s+/)
            if( cols.size() == 2 ) {
                result.put( cols[0], STATUS_MAP.get(cols[1]) )
            }
        }

        return result
    }
}
