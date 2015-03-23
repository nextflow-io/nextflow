/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
/**
 * Executor for ClusterK Cirrus scheduler
 *
 * See http://clusterk.com
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@ServiceName('cirrus')
class CirrusExecutor extends AbstractGridExecutor {

    /**
     * Create a a queue holder for this executor
     * @return
     */
    def TaskMonitor createTaskMonitor() {
        return TaskPollingMonitor.create(session, name, 50, Duration.of('10 sec'))
    }

    def GridTaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir

        final folder = task.workDir
        log.debug "Launching process > ${task.name} -- work folder: $folder"

        final bash = new CirrusWrapperBuilder(task)

        // staging/unstage input/output files
        bash.stagingScript = stagingInputFilesScript(task)
        bash.unstagingScript = unstageOutputFilesScript(task)

        // create the wrapper script
        bash.build()

        return new CirrusTaskHandler(task, this)
    }

    @Override
    String stageInputFileScript( Path path, String targetName ) {
        def op = "es3 -q -v0 --no-stats sync s3:/${path} ."
        if( path.name != targetName )
            op += " && mv ${path.name} ${targetName}"

        return op
    }

    @Override
    String unstageOutputFilesScript( final TaskRun task, final String separatorChar = '\n' ) {

        // collect all the expected names (pattern) for files to be un-staged
        def result = []
        def fileOutNames = task.getOutputFilesNames()
        def normalized = normalizeGlobStarPaths(fileOutNames)

        // create a bash script that will copy the out file to the working directory
        log.trace "Unstaging file path: $normalized"
        if( normalized ) {
            result << ""
            normalized.each {
                result << "es3 -q -v0 --no-stats sync $it s3:/${task.getTargetDir()} || true" // <-- add true to avoid it stops on errors
            }
        }

        return result.join(separatorChar)
    }

    @Override
    protected String getHeaderToken() {
        throw new UnsupportedOperationException()
    }

    @Memoized
    protected List<String> getAwsCredentials() {
        Global.getAwsCredentials()
    }

    @Override
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        // set a tag for the task name
        result << '--tag' << "NAME=${getJobNameFor(task)}"

        result << '--tag' << "uuid=${task.getHashLog()}"

        // Do not propagate host environment to worker nodes
        result << '--no-env'

        def aws = getAwsCredentials()
        if( aws ) {
            result << '-v' << "AWS_ACCESS_KEY_ID=${aws[0]}"
            result << '-v' << "AWS_SECRET_ACCESS_KEY=${aws[1]}"
        }

        // the requested queue name
        if( task.config.queue ) {
            result << '-q' << task.config.queue.toString()
        }

        if( task.config.cpus > 1 ) {
            result << "-c" << task.config.cpus.toString()
        }

        if( task.config.getMemory() ) {
            // convert to MB
            result << '-m' << task.config.getMemory().toMega().toString()
        }

        // -- at the end append the command script wrapped file name
        result.addAll( task.config.getClusterOptionsAsList() )

        return result
    }

    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile) {
        def result = ['ksub']
        getDirectives(task, result)
        result << '-w' << "es3 cat s3:/$scriptFile > ${scriptFile.name} && bash ${scriptFile.name}"
        return result
    }

    @Override
    def parseJobId(String text) {
        return text.trim()
    }

    @Override
    protected List<String> killTaskCommand( jobId ) {
        assert jobId != null
        ['kancel', jobId.toString()]
    }

    @Override
    protected List<String> queueStatusCommand( queue ) {
        def cmd = ['kqueue']
        if( queue )
            cmd << '-q' << queue.toString()

        return cmd
    }

    /*
     *  $ kqueue 3
     *  Task    Attempt Node    Code    State   ExpKb   PeakKb  MaxKb   CPUs    Queue   CmdLine Tags    Priority
     *  3        3      3       -1      DONE    102400  102400  204800  0       default ["/bin/ls"]       0
     */
    @Override
    protected Map<?, QueueStatus> parseQueueStatus(String text) {
        def result = [:]
        text.eachLine { String row, index ->
            if( index< 1 ) return
            def cols = row.trim().split(/\t+/)
            result.put( cols[0], DECODE_STATUS[cols[4]] )
        }
        return result
    }

    private static Map DECODE_STATUS = [
            'PENDING': QueueStatus.PENDING,     // Used only for a Task with dependencies
            'WAITING': QueueStatus.PENDING,     // Task is waiting in queue to be dispatched
            'DISPATCHED': QueueStatus.PENDING,  // Task was dispatched to the node
            'RUNNING': QueueStatus.RUNNING,     // Task is running now
            'CANCELING': QueueStatus.HOLD,      // Task is cancelling now
            'CANCEL_SENT': QueueStatus.HOLD,    // Task cancellation request was
            'DONE': QueueStatus.DONE            // Task completed
    ]
}

@Slf4j
@InheritConstructors
class CirrusTaskHandler extends GridTaskHandler {

    protected ProcessBuilder createProcessBuilder() {

        // -- log the qsub command
        def cli = executor.getSubmitCommandLine(task, wrapperFile)
        log.trace "start process ${task.name} > cli: ${cli}"

        /*
         * launch 'sub' script wrapper
         */
        ProcessBuilder builder = new ProcessBuilder()
                .command( cli as String[] )
                .redirectErrorStream(true)

        return builder
    }

}


/**
 * BASH wrapper builder for running task in Cirrus cluster
 *
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CirrusWrapperBuilder extends BashWrapperBuilder {

    CirrusWrapperBuilder( TaskRun task ) {
        super(task)
        scratch = '$TMPDIR'
        headerScript = fetchScripts(task)
    }

    private String fetchScripts( TaskRun task ) {

        def env = task.workDir.resolve(TaskRun.CMD_ENV)
        def script = task.workDir.resolve(TaskRun.CMD_SCRIPT)
        def infile = task.workDir.resolve(TaskRun.CMD_INFILE)
        def wrapper = task.workDir.resolve(TaskRun.CMD_WRAPPER)

        def ops = []
        ops << '# fetch scripts'
        ops << "es3 test s3:/${env} && es3 cat s3:/${env} > ${env.name}"
        ops << "es3 test s3:/${script} && es3 cat s3:/${script} > ${script.name}"
        ops << "es3 test s3:/${infile} && es3 cat s3:/${infile} > ${infile.name}"
        ops << "es3 test s3:/${wrapper} && es3 cat s3:/${wrapper} > ${wrapper.name}"
        ops << ''

        ops.join('\n')
    }

    @Override
    protected String touchFile( Path file ) {
        "es3 touch s3:/${file}"
    }

    @Override
    protected String fileStr( Path file ) {
        file.getFileName().toString()
    }

    @Override
    protected String copyFile( String name, Path target ) {
        "es3 -q -v 0 --no-stats sync ${name} s3:/${target}"
    }

    protected String exitFile( Path file ) {
        "${file.name} && es3 -q -v 0 --no-stats sync ${file.name} s3:/${file.parent} || true"
    }

}

