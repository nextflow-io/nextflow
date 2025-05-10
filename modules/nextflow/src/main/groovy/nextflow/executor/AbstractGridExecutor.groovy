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

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.processor.TaskConfig
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.Escape
import nextflow.util.Throttle
import org.apache.commons.lang.StringUtils
/**
 * Generic task processor executing a task through a grid facility
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class AbstractGridExecutor extends Executor {

    protected Duration queueInterval

    private final static List<String> INVALID_NAME_CHARS = [ " ", "/", ":", "@", "*", "?", "\\n", "\\t", "\\r", "=" ]

    private Map lastQueueStatus

    /**
     * Initialize the executor class
     */
    @Override
    protected void register () {
        super.register()
        queueInterval = session.getQueueStatInterval(name)
        log.debug "Creating executor '$name' > queue-stat-interval: ${queueInterval}"
    }

    /**
     * Create a a queue holder for this executor
     * @return
     */
    TaskMonitor createTaskMonitor() {
        return TaskPollingMonitor.create(session, name, 100, Duration.of('5 sec'))
    }

    /*
     * Prepare and launch the task in the underlying execution platform
     */
    GridTaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir

        new GridTaskHandler(task, this)
    }

    protected BashWrapperBuilder createBashWrapperBuilder(TaskRun task) {
        // creates the wrapper script
        final builder = new BashWrapperBuilder(task)
        // job directives headers
        builder.headerScript = getHeaderScript(task)
        return builder
    }

    protected String getHeaderScript(TaskRun task) {
        def result = getHeaders(task)
        result += "NXF_CHDIR=${Escape.path(task.workDir)}\n"
        return result
    }

    /**
     * Defines the jobs directive headers
     *
     * @param task
     * @return A multi-line string containing the job directives
     */
    String getHeaders( TaskRun task ) {

        final token = getHeaderToken()
        def result = new StringBuilder()
        def header = new ArrayList(2)
        def dir = getDirectives(task)
        def len = dir.size()-1
        for( int i=0; i<len; i+=2) {
            def opt = dir[i]
            def val = dir[i+1]
            if( opt ) header.add(opt)
            if( val ) header.add(wrapHeader(val))

            if( header ) {
                result << token << ' ' << header.join(' ') << '\n'
            }

            header.clear()
        }

        return result.toString()
    }

    /**
     * @return String used to declare job directives in the job script wrapper
     */
    abstract protected String getHeaderToken()

    /**
     * @param task The current task object
     * @return A list of directives for this task used for the job submission
     */
    final List<String> getDirectives(TaskRun task) {
        getDirectives(task, new ArrayList<String>())
    }

    /**
     * @param task The current task object
     * @param initial An initial list of directives
     * @return A list of directives for this task used for the job submission
     */
    abstract protected List<String> getDirectives(TaskRun task, List<String> initial)

    protected void addClusterOptionsDirective(TaskConfig config, List<String> result) {
        final opts = config.getClusterOptions()
        if( opts instanceof Collection ) {
            for( String it : opts ) {
                result.add(it)
                result.add('')
            }
        }
        else if( opts instanceof CharSequence ) {
            result.add(opts.toString())
            result.add('')
        }
        else if( opts != null ) {
            throw new IllegalArgumentException("Unexpected value for clusterOptions process directive - offending value: $opts")
        }
    }

    /**
     * Given a task returns a *clean* name used to submit the job to the grid engine.
     * That string must not contain blank or special shell characters e.g. parenthesis, etc
     *
     * @param task A {@code TaskRun} instance
     * @return A string that represent to submit job name
     */
    protected String getJobNameFor(TaskRun task) {

        // -- check for a custom `jobName` defined in the nextflow config file
        def customName = resolveCustomJobName(task)
        if( customName )
            return sanitizeJobName(customName)

        // -- if not available fallback on the custom naming strategy

        final result = new StringBuilder("nf-")
        final name = task.getName()
        for( int i=0; i<name.size(); i++ ) {
            final ch = name[i]
            result.append( INVALID_NAME_CHARS.contains(ch) ? "_" : ch )
        }
        return sanitizeJobName(result.toString())
    }

    protected String sanitizeJobName(String name) {
        name.size() > 256 ? name.substring(0,256) : name
    }

    /**
     * Resolve the `jobName` property defined in the nextflow config file
     *
     * @param task
     * @return
     */
    @PackageScope
    String resolveCustomJobName(TaskRun task) {
        try {
            def custom = (Closure)session?.getExecConfigProp(name, 'jobName', null)
            if( !custom )
                return null

            def ctx = [ (TaskProcessor.TASK_CONTEXT_PROPERTY_NAME): task.config ]
            custom.cloneWith(ctx).call()?.toString()
        }
        catch( Exception e ) {
            log.debug "Unable to resolve job custom name", e
            return null
        }
    }

    /**
     * Build up the platform native command line used to submit the job wrapper
     * execution request to the underlying grid, e.g. {@code qsub -q something script.job}
     *
     * @param task The task instance descriptor
     * @return A list holding the command line
     */
    abstract List<String> getSubmitCommandLine(TaskRun task, Path scriptFile)

    /**
     * Defines how script is run the by the grid-engine.
     * @return When {@code true} the launcher script is piped over the submit tool stdin stream,
     *  if {@code false} is specified as an argument on the command line
     */
    protected boolean pipeLauncherScript() { false }

    /**
     * Given the string returned the by grid submit command, extract the process handle i.e. the grid jobId
     */
    abstract parseJobId( String text );

    /**
     * Kill a grid job
     *
     * @param jobId The ID of the job to kill
     */
    void killTask( def jobId )  {
        def cmd = killTaskCommand(jobId)
        def proc = new ProcessBuilder(cmd).redirectErrorStream(true).start()
        proc.waitForOrKill(10_000)
        def ret = proc.exitValue()
        if( ret==0 )
            return

        def m = """\
                Unable to kill pending jobs
                - cmd executed: ${cmd.join(' ')}
                - exit status : $ret
                - output      :
                """
                .stripIndent(true)
        m += proc.text.indent('  ')
        log.debug(m)
    }

    /**
     * The command to be used to kill a grid job
     *
     * @param jobId The job ID to be kill
     * @return The command line to be used to kill the specified job
     */
    protected List<String> killTaskCommand(def jobId) {
        final result = getKillCommand()
        if( jobId instanceof Collection ) {
            jobId.each { result.add(it.toString()) }
            log.trace "Kill command: ${result}"
        }
        else {
            result.add(jobId.toString())
        }
        return result
    }

    /**
     * The command to be used to kill a grid job
     *
     * @return The command line to be used to kill the specified job
     */
    protected abstract List<String> getKillCommand()

    /**
     * Status as returned by the grid engine
     */
    static protected enum QueueStatus { PENDING, RUNNING, HOLD, ERROR, DONE, UNKNOWN }

    /**
     * @return The status for all the scheduled and running jobs
     */
    protected Map<String,QueueStatus> getQueueStatus0(queue) {

        List cmd = queueStatusCommand(queue)
        if( !cmd ) return null

        try {
            log.trace "[${name.toUpperCase()}] getting queue ${queue?"($queue) ":''}status > cmd: ${cmd.join(' ')}"

            final buf = new StringBuilder()
            final process = new ProcessBuilder(cmd).redirectErrorStream(true).start()
            final consumer = process.consumeProcessOutputStream(buf)
            process.waitForOrKill(60_000)
            final exit = process.exitValue(); consumer.join() // <-- make sure sync with the output consume #1045
            final result = buf.toString()

            if( exit == 0 ) {
                log.trace "[${name.toUpperCase()}] queue ${queue?"($queue) ":''}status > cmd exit: $exit\n$result"
                return parseQueueStatus(result)
            }
            else {
                def m = """\
                [${name.toUpperCase()}] queue ${queue?"($queue) ":''}status cannot be fetched
                - cmd executed: ${cmd.join(' ')}
                - exit status : $exit
                - output      :
                """.stripIndent(true)
                m += result.indent('  ')
                log.warn1(m, firstOnly: true)
                return null
            }

        }
        catch( Exception e ) {
            log.warn "[${name.toUpperCase()}] failed to retrieve queue ${queue?"($queue) ":''}status -- See the log file for details", e
            return null
        }
    }

    Map<String,QueueStatus> getQueueStatus(queue) {
        final global = session.getExecConfigProp(name, 'queueGlobalStatus',false)
        if( global ) {
            log.debug1("Executor '$name' fetching queue global status")
            queue = null
        }
        Map<String,QueueStatus> status = Throttle.cache("${name}_${queue}", queueInterval) {
            final result = getQueueStatus0(queue)
            log.trace "[${name.toUpperCase()}] queue ${queue?"($queue) ":''}status >\n" + dumpQueueStatus(result)
            return result
        }
        // track the last status for debugging purpose
        lastQueueStatus = status
        return status
    }

    @PackageScope
    final String dumpQueueStatus(Map<String,QueueStatus> statusMap) {
        if( statusMap == null )
            return '  (null)'
        if( statusMap.isEmpty() )
            return '  (empty)'

        def result = new StringBuilder()
        statusMap?.each { k, v ->
            result << '  job: ' << StringUtils.leftPad(k?.toString(),6) << ': ' << v?.toString() << '\n'
        }
        return result.toString()
    }

    @PackageScope
    final String dumpQueueStatus() {
        dumpQueueStatus(lastQueueStatus)
    }

    /**
     * @param queue The command for which the status of jobs has be to read
     * @return The command line to be used to retried the job statuses
     */
    protected abstract List<String> queueStatusCommand(queue)

    /**
     * Parse the command stdout produced by the command line {@code #queueStatusCommand}
     * @param text
     * @return
     */
    protected abstract Map<String,QueueStatus> parseQueueStatus( String text )

    boolean checkStartedStatus(jobId, queueName) {
        assert jobId

        // -- fetch the queue status
        final queue = getQueueStatus(queueName)
        if( !queue )
            return false

        if( !queue.containsKey(jobId) )
            return false
        if( queue.get(jobId) == QueueStatus.PENDING )
            return false

        return true
    }

    /**
     * Verify that a job in a 'active' state i.e. RUNNING or HOLD
     *
     * @param jobId The job for which verify the status
     * @return {@code true} if the job is in RUNNING or HOLD status, or even if it is temporarily unable
     *  to retrieve the job status for some
     */
    boolean checkActiveStatus( jobId, queue ) {
        assert jobId

        // -- fetch the queue status
        final status = getQueueStatus(queue)

        if( status == null ) { // no data is returned, so return true
            return true
        }

        if( !status.containsKey(jobId) ) {
            log.trace "[${name.toUpperCase()}] queue ${queue?"($queue) ":''}status > map does not contain jobId: `$jobId`"
            return false
        }

        final result = status[jobId.toString()] == QueueStatus.RUNNING || status[jobId.toString()] == QueueStatus.HOLD
        log.trace "[${name.toUpperCase()}] queue ${queue?"($queue) ":''}status > jobId `$jobId` active status: $result"
        return result
    }

    protected String wrapHeader( String str ) { str }

    protected String quote(Path path) {
        def str = Escape.path(path)
        path.toString() != str ? "\"$str\"" : str
    }

    @Override
    boolean isContainerNative() {
        // when fusion is enabled it behaves as a native container environment
        // because the command wrapper script should not manage the container execution.
        // Instead, it is the command wrapper script that is launched run within a container process.
        return isFusionEnabled()
    }
}

