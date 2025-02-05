/*
 * Copyright 2013-2024, Seqera Labs
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
import java.nio.file.Paths
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.fusion.FusionHelper
import nextflow.processor.TaskArrayRun
import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
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
@CompileStatic
class LsfExecutor extends AbstractGridExecutor implements TaskArrayExecutor {

    static private Pattern KEY_REGEX = ~/^[A-Z_0-9]+=.*/

    static private Pattern QUOTED_STRING_REGEX = ~/"((?:[^"\\]|\\.)*)"(\s*#.*)?/

    private boolean perJobMemLimit

    private boolean perTaskReserve

    /*
     * If LSF_UNIT_FOR_LIMITS is not defined in lsf.conf, then the default setting is in KB, and for RUSAGE it is MB
     * see https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.3/lsf_config_ref/lsf.conf.lsf_unit_for_limits.5.html
     */
    private String memUnit = 'KB'

    private String usageUnit = 'MB'

    protected boolean getPerJobMemLimit() { perJobMemLimit }
    protected boolean getPerTaskReserve() { perTaskReserve }
    protected String getMemUnit() { memUnit }
    protected String getUsageUnit() { usageUnit }
    
    /**
     * Gets the directives to submit the specified task to the cluster for execution
     *
     * @param task A {@link TaskRun} to be submitted
     * @param result The {@link List} instance to which add the job directives
     * @return A {@link List} containing all directive tokens and values.
     */
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        result << '-o' << (task.isArray() ?  '/dev/null' : task.workDir.resolve(TaskRun.CMD_LOG).toString())

        // add other parameters (if any)
        if( task.config.queue ) {
            result << '-q'  << (task.config.queue as String)
        }

        //number of cpus for multiprocessing/multi-threading
        if( task.config.getCpus() > 1 ) {
            result << "-n" << task.config.getCpus().toString()
            result << "-R" << "span[hosts=1]"
        }

        if( task.config.getTime() ) {
            result << '-W' << task.config.getTime().format('HH:mm')
        }

        if( task.config.getMemory() ) {
            def mem = task.config.getMemory()
            // LSF mem limit can be both per-process and per-job
            // depending a system configuration setting -- see https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.3/lsf_config_ref/lsf.conf.lsb_job_memlimit.5.dita
            // When per-process is used (default) the amount of requested memory
            // is divided by the number of used cpus (processes)
            def mem1 = ( task.config.getCpus() > 1 && !perJobMemLimit ) ? mem.div(task.config.getCpus() as int) : mem
            def mem2 = ( task.config.getCpus() > 1 && perTaskReserve ) ? mem.div(task.config.getCpus() as int) : mem

            result << '-M' << String.valueOf(mem1.toUnit(memUnit))
            result << '-R' << "select[mem>=${mem.toUnit(memUnit)}] rusage[mem=${mem2.toUnit(usageUnit)}]".toString()
        }

        def disk = task.config.getDisk()
        if( disk ) {
            result << '-R' << "select[tmp>=${disk.toUnit(memUnit)}] rusage[tmp=${disk.toUnit(usageUnit)}]".toString()
        }

        // -- the job name
        if( task instanceof TaskArrayRun ) {
            final arraySize = task.getArraySize()
            result << '-J' << "${getJobNameFor(task)}[1-${arraySize}]".toString()
        }
        else {
            result << '-J' << getJobNameFor(task)
        }

        // -- at the end append the command script wrapped file name
        addClusterOptionsDirective(task.config, result)

        // add account from config
        final account = session.getExecConfigProp(getName(), 'account', null) as String
        if( account ) {
            result << '-G' << account
        }

        return result
    }

    @Override
    protected void addClusterOptionsDirective(TaskConfig config, List<String> result) {
        final opts = config.getClusterOptions()
        // when cluster options are defined as a list rely on default behavior
        if( opts instanceof Collection ) {
            super.addClusterOptionsDirective(config,result)
        }
        // when cluster options are a string value use the `getClusterOptionsAsList` for backward compatibility
        else if( opts instanceof CharSequence ) {
            result.addAll( config.getClusterOptionsAsList() )
        }
        else if( opts != null ) {
            throw new IllegalArgumentException("Unexpected value for clusterOptions process directive - offending value: $opts")
        }
    }

    @Override
    String sanitizeJobName( String name ) {
        // LSF does not allow square brackets in job names except for job arrays
        name = name.replace('[','').replace(']','')
        // Old LSF versions do not allow job names longer than 511 chars
        name.size()>511 ? name.substring(0,511) : name
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
        // note: use the `-w` option to avoid that the printed jobid may be truncated when exceed 7 digits
        // see https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.3/lsf_config_ref/lsf.conf.lsb_jobid_disp_length.5.html
        final result = ['bjobs', '-w']

        if( queue )
            result << '-q' << queue.toString()

        return result
    }

    // https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.3/lsf_command_ref/bjobs.zz4category.description.1.html
    private static Map<String,QueueStatus> DECODE_STATUS = [
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

        final Map<String, QueueStatus> result = new LinkedHashMap<>()
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

            final jobId = line.tokenize(' ')[col1]
            final status = line.tokenize(' ')[col2]
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

    protected String getEnv0(String name) {
        System.getenv(name)
    }

    /**
     * parse the memory units from $LSF_ENVDIR/lsf.conf,
     * which should be accessible via all nodes. Otherwise, default is KB
     * -- see https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.3/lsf_config_ref/lsf.conf.lsf_unit_for_limits.5.html
     */
    protected Map<String,String> parseLsfConfig() {
        final Map<String,String> result = new LinkedHashMap<>(20)

        // check environment variable exists
        def envDir = getEnv0('LSF_ENVDIR')
        if( !envDir )
            return result

        def envFile = Paths.get(envDir).resolve("lsf.conf")
        if( !envFile.exists() )
            return result

        for( def line : envFile.readLines() ) {
            if( !KEY_REGEX.matcher(line).matches() )
                continue
            def entry = line.tokenize('=')
            if( entry.size() != 2 )
                continue
            def key = entry[0]
            def value = entry[1]
            def matcher = QUOTED_STRING_REGEX.matcher(value)
            if( matcher.matches() ) {
                value = matcher.group(1)
            }
            else {
                int p = value.indexOf('#')
                value = p==-1 ? value.trim() : value.substring(0,p).trim()
            }

            result.putAt(key,value)
        }

        return result
    }

    @Override
    void register() {
        super.register()
        final conf = parseLsfConfig()

        // lsf mem unit
        // https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.3/lsf_config_ref/lsf.conf.lsf_unit_for_limits.5.html
        if( conf.get('LSF_UNIT_FOR_LIMITS') ) {
            memUnit = usageUnit = conf.get('LSF_UNIT_FOR_LIMITS')
            log.debug "[LSF] Detected lsf.conf LSF_UNIT_FOR_LIMITS=$memUnit"
        }
        
        // per job mem limit
        // https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.3/lsf_config_ref/lsf.conf.lsb_job_memlimit.5.dita
        if( conf.get('LSB_JOB_MEMLIMIT') ) {
            final str = conf.get('LSB_JOB_MEMLIMIT').toUpperCase()
            perJobMemLimit = str == 'Y'
            log.debug "[LSF] Detected lsf.conf LSB_JOB_MEMLIMIT=$str ($perJobMemLimit)"
        }

        perJobMemLimit = session.getExecConfigProp(name, 'perJobMemLimit', perJobMemLimit)

        // per task reserve https://github.com/nextflow-io/nextflow/issues/1071#issuecomment-481412239
        if( conf.get('RESOURCE_RESERVE_PER_TASK') ) {
            final str = conf.get('RESOURCE_RESERVE_PER_TASK').toUpperCase()
            perTaskReserve = str == 'Y'
            log.debug "[LSF] Detected lsf.conf RESOURCE_RESERVE_PER_TASK=$str ($perTaskReserve)"
        }

        perTaskReserve = session.getExecConfigProp(name, 'perTaskReserve', perTaskReserve)
    }

    @Override
    boolean isFusionEnabled() {
        return FusionHelper.isFusionEnabled(session)
    }

    @Override
    String getArrayIndexName() {
        return 'LSB_JOBINDEX'
    }

    @Override
    int getArrayIndexStart() {
        return 1
    }

    @Override
    String getArrayTaskId(String jobId, int index) {
        assert jobId, "Missing 'jobId' argument"
        return "${jobId}[${index + 1}]"
    }

}
