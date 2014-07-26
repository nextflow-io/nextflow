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
import static nextflow.processor.TaskHandler.Status.COMPLETED
import static nextflow.processor.TaskHandler.Status.RUNNING
import static nextflow.processor.TaskHandler.Status.SUBMITTED

import java.nio.file.Path

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessFailedException
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.trace.TraceRecord
import nextflow.util.Duration
import org.ggf.drmaa.ExitTimeoutException
import org.ggf.drmaa.JobInfo
import org.ggf.drmaa.JobTemplate
import org.ggf.drmaa.Session as DrmaaSession
import org.ggf.drmaa.SessionFactory

/**
 * A job executor based on DRMAA interface
 *
 * @link http://www.drmaa.org
 * @link http://gridscheduler.sourceforge.net/howto/drmaa_java.html
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class DrmaaExecutor extends Executor {

    static private DrmaaSession drmaa

    /**
     * Create a DRMAA session object the very first time an executor is created
     */
    void register() {
        drmaa = SessionFactory.getFactory().getSession()
        drmaa.init('')
        session.onShutdown { drmaa.exit() }
    }

    /**
     * @return The DRMAA session instance
     */
    DrmaaSession getDrmaaSession() { drmaa }

    /**
     * Create a a queue holder for this executor
     * @return
     */
    def TaskMonitor createTaskMonitor() {
        return TaskPollingMonitor.create(session, name, 50, Duration.of('1 sec'))
    }

    /*
     * Prepare and launch the task in the underlying execution platform
     */
    def DrmaaTaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir

        final folder = task.workDir
        log.debug "Launching process > ${task.name} -- work folder: $folder"

        final bash = new BashWrapperBuilder(task)

        // staging/unstage input/output files
        bash.stagingScript = stagingInputFilesScript(task)
        bash.unstagingScript = unstageOutputFilesScript(task)

        // create the wrapper script
        bash.build()

        return new DrmaaTaskHandler(task, taskConfig, this)
    }


}


/**
 * Task handler for DRMAA jobs
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class DrmaaTaskHandler extends TaskHandler {

    /** The target executor platform */
    final DrmaaExecutor executor

    /** Location of the file created when the job is started */
    final Path startFile

    /** Location of the file created when the job is terminated */
    final Path exitFile

    /** Location of the file holding the task std output */
    final Path outputFile

    /** The wrapper file used to execute the user script */
    final Path wrapperFile

    /** The unique job ID as provided by the underlying grid platform */
    private String jobId

    private String taskName

    private File workDir

    private DrmaaSession drmaa

    private JobInfo fJobInfo

    DrmaaTaskHandler(TaskRun task, TaskConfig taskConfig, DrmaaExecutor executor) {
        super(task, taskConfig)

        this.executor = executor
        this.drmaa = executor.getDrmaaSession()
        this.startFile = task.workDir.resolve(TaskRun.CMD_START)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.wrapperFile = task.workDir.resolve(TaskRun.CMD_RUN)
        this.taskName = "nf-${task.name.replace(' ','_')}"
        this.workDir = task.workDir?.toFile()
    }

    @PackageScope
    String getOptions() {

        List result = []
        result << '-notify'

        // the requested queue name
        if( taskConfig.queue ) {
            result << '-q'  << taskConfig.queue
        }

        // max task duration
        if( taskConfig.maxDuration ) {
            final duration = taskConfig.maxDuration as Duration
            result << "-l" << "h_rt=${duration.format('HH:mm:ss')}"
        }

        // task max memory
        if( taskConfig.maxMemory ) {
            result << "-l" << "virtual_free=${taskConfig.maxMemory.toString().replaceAll(/[\sB]/,'')}"
        }

        // -- at the end append the command script wrapped file name
        if( taskConfig.clusterOptions ) {
            result << taskConfig.clusterOptions
        }

        result << '-b y'
        return result.join(' ')
    }


    @Override
    void submit() {

        JobTemplate template = drmaa.createJobTemplate()
        template.setWorkingDirectory(workDir.toString())
        template.setRemoteCommand('/bin/bash')
        template.setArgs( [wrapperFile.toString()] )
        template.setJoinFiles(true)
        template.setOutputPath( ':/dev/null' )
        template.setNativeSpecification( getOptions() )

        try {
            jobId = drmaa.runJob(template)
            this.status = SUBMITTED
        }
        catch( Exception e ) {
            status = COMPLETED
            throw new ProcessFailedException("DRMAA exception submitting process '${task.name}' for execution", e )
        }
        finally {
            drmaa.deleteJobTemplate(template)
        }

    }

    @Override
    boolean checkIfRunning() {

        if( isSubmitted()  ) {

            if( startFile && startFile.exists() && startFile.lastModified() > 0) {
                status = RUNNING
                return true
            }

        }

        return false
    }


    private getJobInfo() {
        if( fJobInfo )
            return fJobInfo

        try {
            return fJobInfo = drmaa.wait(jobId, 0)
        }
        catch( ExitTimeoutException e ) {
            return null
        }
    }

    @Override
    boolean checkIfCompleted() {

        def job
        if( isRunning() && (job=getJobInfo()) != null) {

            if( job.wasAborted() || job.hasSignaled() ) {
                status = COMPLETED
                return true
            }
            else if ( job.hasExited() && exitFile.lastModified() > 0 ) {
                status = COMPLETED
                task.stdout = outputFile
                task.exitStatus = job.getExitStatus()
                return true
            }

        }

        return false
    }

    @Override
    void kill() {
        drmaa.control(jobId, DrmaaSession.TERMINATE);
    }

    /*

     Example content of the {@link JobInfo#getResourceUsage()} map

     mem = 0.0000
     io = 0.0000
     iow = 0.0000
     maxvmem = 0.0000
     signal = 0.0000
     vmem = 0.0000
     cpu = 0.0040
     priority = 0.0000
     exit_status = 0.0000
     start_time = 1406265009.0000
     submission_time = 1406264935.0000
     end_time = 1406265009.0000
     ru_utime = 0.0020
     ru_maxrss = 5392.0000
     ru_idrss = 0.0000
     ru_inblock = 8.0000
     ru_nivcsw = 0.0000
     ru_isrss = 0.0000
     ru_oublock = 8.0000
     ru_majflt = 0.0000
     ru_minflt = 647.0000
     ru_msgrcv = 0.0000
     ru_nswap = 0.0000
     ru_ixrss = 0.0000
     ru_msgsnd = 0.0000
     ru_wallclock = 0.0000
     ru_nsignals = 0.0000
     ru_stime = 0.0020
     ru_nvcsw = 18.0000
     ru_ismrss = 0.0000
     acct_mem = 0.0000
     acct_iow = 0.0000
     acct_io = 0.0000
     acct_maxvmem = 0.0000
     acct_cpu = 0.0040

     */
    @Override
    public TraceRecord getTraceRecord() {
        Map<String,String> resources = getJobInfo()?.getResourceUsage() ?: [:]

        def trace = super.getTraceRecord()
        trace.nativeId = jobId
        trace.exit = task.exitStatus

        if( resources ) {
            if( resources.containsKey('submission_time'))
                trace.submit = resources['submission_time']?.trimDotZero()?.toLong()

            if( resources.containsKey('start_time'))
                trace.start = resources['start_time']?.trimDotZero()?.toLong()

            if( resources.containsKey('end_time'))
                trace.complete = resources['end_time']?.trimDotZero()?.toLong()

            if( resources.containsKey('maxvmem') )
                trace.mem = resources['maxvmem']

            if( resources.containsKey('cpu'))
                trace.cpu = resources['cpu']

            if( resources.containsKey('io'))
                trace.io = resources['io']
        }

        return trace
    }


}
