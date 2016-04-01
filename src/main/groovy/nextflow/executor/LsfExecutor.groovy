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

import nextflow.processor.TaskRun
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

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        result << '-cwd' << task.workDir.toString()
        result << '-o' << task.workDir.resolve(TaskRun.CMD_LOG).toString()

        // add other parameters (if any)
        if( task.config.queue ) {
            result << '-q'  << (task.config.queue as String)
        }

        //number of cpus for multiprocessing/multi-threading
        if( task.config.cpus > 1 ) {
            result << "-n" << task.config.cpus.toString()
            result << "-R" << "span[hosts=1]"
        }

        if( task.config.time ) {
            result << '-W' << task.config.getTime().format('HH:mm')
        }

        if( task.config.getMemory() ) {
            def mem = task.config.getMemory()
            // LSF specify per-process (per-core) memory limit (in MB)
            if( task.config.cpus > 1 ) {
                long bytes = mem.toBytes().intdiv(task.config.cpus as int)
                mem = new MemoryUnit(bytes)
            }
            // convert to MB
            result << '-M' << String.valueOf(mem.toMega())
            result << '-R' << "rusage[mem=${mem.toMega()}]"
        }

        // -- the job name
        result << '-J' << getJobNameFor(task)

        // -- at the end append the command script wrapped file name
        result.addAll( task.config.getClusterOptionsAsList() )

        return result
    }


    /**
     * The command line to submit this job
     *
     * @param task The {@link TaskRun} instance to submit for execution to the cluster
     * @param scriptFile The file containing the job launcher script
     * @return A list representing the submit command line
     */
    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) { ['bsub'] }

    /**
     * @return {@code true} since BSC grid requires the script to be piped to the {@code bsub} command
     */
    @Override
    protected boolean pipeLauncherScript() { true }

    protected String getHeaderToken() { '#BSUB' }


    /**
     * Parse the string returned by the {@code bsub} command and extract the job ID string
     *
     * @param text The string returned when submitting the job
     * @return The actual job ID string
     */
    @Override
    def parseJobId(String text) {

        def pattern = ~/Job <(\d+)> is submitted/
        for( String line : text.readLines() ) {
            def m = pattern.matcher(line)
            if( m.find() ) {
                return m.group(1)
            }
        }

        new IllegalStateException("Invalid LSF submit response:\n$text\n\n");
    }

    @Override
    String getKillCommand() { 'bkill' }

    @Override
    protected List<String> queueStatusCommand( queue ) {

        def result = ['bjobs', '-o',  'JOBID STAT SUBMIT_TIME delimiter=\',\'', '-noheader']

        if( queue )
            result << '-q' << queue.toString()

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


    private static String SPECIAL_CHARS = ' []|&!<>'

    @Override
    protected String wrapHeader( String str ) {
        boolean needWrap = false

        if( !str ) return str

        for( int i=0; i<SPECIAL_CHARS.size(); i++ ) {
            for( int x=0; x<str.size(); x++ ) {
                if( str.charAt(x) == SPECIAL_CHARS.charAt(i) ) {
                    needWrap = true
                    break
                }
            }
        }

        return needWrap ? "\"$str\"" : str
    }

}

