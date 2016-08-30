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

import com.upplication.s3fs.S3Path
import groovy.transform.InheritConstructors
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.exception.AbortOperationException
import nextflow.processor.TaskBean
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.ServiceName

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

    @Override
    void register() {
        super.register()
        // check that the working directory is on S3 storage
        if( !(session.workDir instanceof S3Path) ) {
            session.abort()
            throw new AbortOperationException("When using `cirrus` executor a S3 bucket must be provided as working directory -- Add the option `-w s3://<your-bucket/path>` to your run command line")
        }
    }

    /**
     * Create a a queue holder for this executor
     * @return
     */
    @Override
    def TaskMonitor createTaskMonitor() {
        return TaskPollingMonitor.create(session, name, 50, Duration.of('10 sec'))
    }

    @Override
    def GridTaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir

        return new CirrusTaskHandler(task, this)
    }

    @Override
    protected BashWrapperBuilder createBashWrapperBuilder(TaskRun task) {
        final folder = task.workDir
        log.debug "Launching process > ${task.name} -- work folder: $folder"

        // create the wrapper script
        final bean = new TaskBean(task)
        final bash = new BashWrapperBuilder(bean, new CirrusFileCopyStrategy(bean))
        bash.scratch = '$TMPDIR'
        return bash
    }

    @Override
    protected String getJobNameFor(TaskRun task) {
        task.getName().replace(' ','_')
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

        // max task duration
        if( task.config.time ) {
            final time = task.config.getTime().toSeconds()
            result << '--timeout' << String.valueOf(time)
        }

        if( task.config.getMemory() ) {
            // convert to MB
            result << '-m' << (task.config.getMemory().toMega() as String)
        }

        if( task.config.getDisk() ) {
            result << '-dm' << (task.config.getDisk().toMega() as String)
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
    protected String getKillCommand() { 'kancel' }

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



