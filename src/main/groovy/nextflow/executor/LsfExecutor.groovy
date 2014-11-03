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

import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.MemoryUnit


/**
 * Processor for LSF resource manager (DRAFT)
 *
 * See
 * http://en.wikipedia.org/wiki/Platform_LSF
 * https://doc.zih.tu-dresden.de/hpc-wiki/bin/view/Compendium/PlatformLSF
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LsfExecutor extends AbstractGridExecutor {

    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) {

        // note: LSF requires the job script file to be executable
        scriptFile.setPermissions(7,0,0)

        final result = new ArrayList<String>()
        result << 'bsub'
        result << '-cwd' << task.workDir?.toString()
        result << '-o' << '/dev/null'

        // add other parameters (if any)
        if( taskConfig.queue ) {
            result << '-q'  << (taskConfig.queue as String)
        }

        //number of cpus for multiprocessing/multi-threading
        if( taskConfig.cpus ) {
            result << "-n" << taskConfig.cpus.toString()
            result << "-R" << "span[hosts=1]"
        }

        if( taskConfig.time ) {
            result << '-W' << (taskConfig.time as Duration).format('HH:mm')
        }

        if( taskConfig.getMemory() ) {
            def mem = taskConfig.getMemory()
            // LSF specify per-process (per-core) memory limit (in MB)
            if( taskConfig.cpus && taskConfig.cpus>1 ) {
                long bytes = mem.toBytes() / (int)taskConfig.cpus
                mem = new MemoryUnit(bytes)
            }
            // convert to MB
            result << '-M' << mem.toMega().toString()
        }

        // -- the job name
        result << '-J' << getJobNameFor(task)

        // -- at the end append the command script wrapped file name
        if( taskConfig.clusterOptions ) {
            result.addAll( getClusterOptionsAsList() )
        }

        // -- last entry to 'script' file name
        result << "./${scriptFile.getName()}"

        return result

    }

    @Override
    def parseJobId(String text) {

        def pattern = ~/Job <(\d+)> is submitted/
        for( String line : text.readLines() ) {
            def m = pattern.matcher(line)
            if( m.find() ) {
                return m[0][1].toString()
            }
        }

        new IllegalStateException("Invalid LSF submit response:\n$text\n\n");
    }

    @Override
    List<String> killTaskCommand( def jobId ) {
        ['bkill', jobId?.toString() ]
    }

    @Override
    protected List<String> queueStatusCommand(Object queue) {

        def result = ['bjobs', '-o',  'JOBID STAT SUBMIT_TIME delimiter=\',\'', '-noheader']

        if( queue )
            result << '-q' << queue

        return result

    }

    private static Map DECODE_STATUS = [
            'PEND': QueueStatus.PENDING,
            'RUN': QueueStatus.RUNNING,
            'PSUSP': QueueStatus.HOLD,
            'USUSP': QueueStatus.HOLD,
            'SSUSP': QueueStatus.HOLD,
            'DONE': QueueStatus.DONE,
            'EXIT': QueueStatus.ERROR,
            'UNKWN': QueueStatus.ERROR,
            'ZOMBI': QueueStatus.ERROR,
    ]

    @Override
    protected Map<?, QueueStatus> parseQueueStatus(String text) {

        def result = [:]

        text.eachLine { String line ->
            def cols = line.split(',')
            if( cols.size() == 3 ) {
                result.put( cols[0], DECODE_STATUS.get(cols[1]) )
            }
        }

        return result
    }
}

