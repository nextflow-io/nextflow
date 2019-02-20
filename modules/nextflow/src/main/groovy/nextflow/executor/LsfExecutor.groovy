/*
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

import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
import nextflow.util.MemoryUnit
/**
 * Processor for LSF resource manager
 *
 * See
 * http://en.wikipedia.org/wiki/Platform_LSF
 * https://doc.zih.tu-dresden.de/hpc-wiki/bin/view/Compendium/PlatformLSF
 * https://www.ibm.com/support/knowledgecenter/SSETD4/product_welcome_lsf.html
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class LsfExecutor extends AbstractGridExecutor {

    private static char BLANK = ' ' as char

    /**
    * parse the memory units from $LSF_ENVDIR/lsf.conf, 
    * which should be accessible via all nodes. Otherwise, default is KB
    * -- see https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.3/lsf_config_ref/lsf.conf.lsf_unit_for_limits.5.html
    */
    def getLSFmemoryUnits() {
        // check environment variable exists
        def envDir = System.getenv("LSF_ENVDIR") 
        if ( !envDir ) return null
        def envFile = new File(envDir, "lsf.conf")
        if ( !envFile.exists() ) return null
        for (def line : envFile.readLines() ){
            if ( line.startsWith("LSF_UNIT_FOR_LIMITS=") ){
                def memoryUnits = line.split("LSF_UNIT_FOR_LIMITS=")[1]
            }
        }
        return(memoryUnits)
    }


    /**
     * Add an init method in the LsfExecutor overriding the base init() in AbstractGridExecutor.groovy
    */
    @Override
    void init() {
        super.init()
        queueInterval = session.getQueueStatInterval(name)
        log.debug "Creating executor '$name' > queue-stat-interval: ${queueInterval}"
        def memoryUnitLSF = getLSFmemoryUnits()
    }

    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {

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

            if ( memoryUnitLSF ){   // if LSF_UNIT_FOR_LIMITS is set by $LSF_ENVDIR/lsf.conf and accessible
                // parse this unit and use this
                def mem = task.config.getMemory()
                // LSF mem limit can be both per-process and per-job
                if( task.config.cpus > 1 && !perJobMemLimit ) {
                    long bytes = mem.toBytes().intdiv(task.config.cpus as int)
                    result << '-M' << String.valueOf(MemoryUnit.of(bytes).toUnit(memoryUnitLSF))
                }
                else {
                    result << '-M' << String.valueOf(mem.toUnit(memoryUnitLSF))
                }
                result << '-R' << "select[mem>=${mem.toUnit(memoryUnitLSF)}] rusage[mem=${mem.toUnit(memoryUnitLSF)}]".toString()
            } else{
                // use default KB
                def mem = task.config.getMemory()
                // LSF mem limit can be both per-process and per-job
                if( task.config.cpus > 1 && !perJobMemLimit ) {
                    long bytes = mem.toBytes().intdiv(task.config.cpus as int)
                    result << '-M' << String.valueOf(MemoryUnit.of(bytes).toKilo())
                }
                else {
                    result << '-M' << String.valueOf(mem.toKilo())
                }
                result << '-R' << "select[mem>=${mem.toKilo()}] rusage[mem=${mem.toKilo()}]".toString()
            }
        }

        def disk = task.config.getDisk()
        if( disk ) {
            if ( memoryUnitLSF ){ 
                result << '-R' << "select[tmp>=$disk.toUnit(memoryUnitLSF)] rusage[tmp=$disk.toUnit(memoryUnitLSF)]".toString()
            }
            else {
                result << '-R' << "select[tmp>=$disk.kilo] rusage[tmp=$disk.kilo]".toString()   // using the default as kilo              
            }
        }

        // -- the job name
        result << '-J' << getJobNameFor(task)

        // -- at the end append the command script wrapped file name
        result.addAll( task.config.getClusterOptionsAsList() )

        return result
    }

    protected boolean isPerJobMemLimit() {
        session.getExecConfigProp(name, 'perJobMemLimit', false)
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

        new IllegalStateException("[LSF] Invalid submit response:\n$text\n\n");
    }

    @Override
    protected List<String> getKillCommand() { ['bkill'] }

    @Override
    protected List<String> queueStatusCommand( queue ) {
        // note: use the `-w` option to avoid that the printed jobid maybe truncated when exceed 7 digits
        // see https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.3/lsf_config_ref/lsf.conf.lsb_jobid_disp_length.5.html
        final result = ['bjobs', '-w']

        if( queue )
            result << '-q' << queue.toString()

        return result
    }

    // https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.3/lsf_command_ref/bjobs.zz4category.description.1.html
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
    protected Map<String, QueueStatus> parseQueueStatus(String text) {

        def result = [:]
        def col1 = -1
        def col2 = -1
        for( String line : text.readLines() ) {
            if( !line )
                continue

            if( col1==-1 && col2==-1 ) {
                col1 = line.tokenize(' ').indexOf('JOBID')
                col2 = line.tokenize(' ').indexOf('STAT')
                continue
            }

            def jobId = line.tokenize(' ')[col1]
            def status = line.tokenize(' ')[col2]
            if( jobId )
                result[jobId] = DECODE_STATUS.get(status)
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
