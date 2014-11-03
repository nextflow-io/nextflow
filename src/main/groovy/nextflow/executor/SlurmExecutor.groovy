/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
import nextflow.util.Duration

/**
 * Processor for SLURM resource manager (DRAFT)
 *
 * See http://computing.llnl.gov/linux/slurm/
 *
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@InheritConstructors
class SlurmExecutor extends AbstractGridExecutor {


    /**
     * -c, --cpus-per-task=ncpus   number of cpus required per task
     * -D, --chdir=path            change remote current working directory
     * -e, --error=err             location of stderr redirection
     * -E, --preserve-env          env vars for node and task counts override
     * -i, --input=in              location of stdin redirection
     * -J, --job-name=jobname      name of job
     * -o, --output=out            location of stdout redirection
     * -Q, --quiet                 quiet mode (suppress informational messages)
     * -t, --time=minutes          time limit
     * -u, --unbuffered            do not line-buffer stdout/err
     * --mem=MB                minimum amount of real memory
     * --mincpus=n             minimum number of logical processors (threads) per node
     * --tmp=MB                minimum amount of temporary disk
     * --mem-per-cpu=MB        maximum amount of real memory per allocated cpu required by the job. --mem >= --mem-per-cpu if --mem is specified.
     *
     * @param task
     * @return
     */
    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {

        final result = new ArrayList<String>()

        result << 'sbatch'
        result << '-D' << task.workDir.toString()
        result << '-J' << getJobNameFor(task)
        result << '-o' << '/dev/null'

        if( taskConfig.cpus ) {
            result << '-c' << taskConfig.cpus.toString()
        }

        if( taskConfig.time ) {
            result << '-t' << (taskConfig.time as Duration).format('HH:mm:ss')
        }

        if( taskConfig.getMemory() ) {
            result << '--mem' << taskConfig.getMemory().toMega().toString()
        }

        // -- at the end append the command script wrapped file name
        if( taskConfig.clusterOptions ) {
            result.addAll( getClusterOptionsAsList() )
        }

        // -- last entry to 'script' file name
        // replace with the 'shell' attribute
        result << scriptFile.getName()

    }

    @Override
    def parseJobId(String text) {
        def pattern = ~ /Submitted batch job (\d+)/
        for( String line : text.readLines() ) {
            def m = pattern.matcher(line)
            if( m.matches() ) {
                return m[0][1].toString()
            }
        }

        throw new IllegalStateException("Invalid SLURM submit response:\n$text\n\n")
    }


    protected List<String> killTaskCommand(def jobId) {
        ['scancel', jobId?.toString() ]
    }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        if( queue )
            log.debug "SLURM executor does not support queue parameter on queue status"

        return ['squeue','-h','-o','%i %t']
    }

    static private Map STATUS_MAP = [
            'PD': QueueStatus.PENDING,  // (pending)
            'R': QueueStatus.RUNNING,   // (running)
            'CA': QueueStatus.ERROR,    // (cancelled)
            'CF': QueueStatus.PENDING,  // (configuring)
            'CG': QueueStatus.RUNNING,  // (completing)
            'CD': QueueStatus.DONE,     // (completed)
            'F': QueueStatus.ERROR,     // (failed),
            'TO': QueueStatus.ERROR,    // (timeout),
            'NF': QueueStatus.ERROR     // (node failure)
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