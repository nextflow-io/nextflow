/*
 * Copyright (c) 2012, the authors.
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

import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
import nextflow.util.ByteDumper
import nextflow.util.Duration
import org.apache.commons.io.IOUtils
import org.ggf.drmaa.JobTemplate
import org.ggf.drmaa.SessionFactory

/**
 * Use the DRMAA Java API to manage the submission to the underlying resource manager
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
class DrmaaExecutor extends AbstractExecutor {

    private static final COMMAND_OUT_FILENAME = '.command.out'
    private static final COMMAND_RUNNER_FILENAME = '.command.run'
    private static final COMMAND_ENV_FILENAME = '.command.env'
    private static final COMMAND_SCRIPT_FILENAME = '.command.sh'
    protected static final JOB_OUT_FILENAME = '.job.out'
    protected static final JOB_SCRIPT_FILENAME = '.job.run'

    @Lazy
    private static factory = SessionFactory.getFactory();

    @Override
    void launchTask(TaskRun task) {
        assert task
        assert task.@script
        assert task.workDirectory

        final scratch = task.workDirectory
        log.debug "Lauching task > ${task.name} -- scratch folder: $scratch"

        /*
         * save the environment to a file
         */
        final envMap = task.processor.getProcessEnvironment()
        final envBuilder = new StringBuilder()
        envMap.each { name, value ->
            if( name ==~ /[a-zA-Z_]+[a-zA-Z0-9_]*/ ) {
                envBuilder << "export $name='$value'" << '\n'
            }
            else {
                log.debug "Task ${task.name} > Invalid environment variable name: '${name}'"
            }
        }
        new File(scratch, COMMAND_ENV_FILENAME).text = envBuilder.toString()

        /*
         * save the main script file
         */
        File cmdOutFile = new File(scratch, COMMAND_OUT_FILENAME)
        def scriptFile = new File(scratch, COMMAND_SCRIPT_FILENAME)
        scriptFile.text = task.processor.normalizeScript(task.script.toString())
        scriptFile.setExecutable(true)

        def runnerText = """
                    source ${COMMAND_ENV_FILENAME}
                    [ ! -z \$TMPDIR ] && cd \$TMPDIR
                    ./${COMMAND_SCRIPT_FILENAME}
                    """
        def runnerFile = new File(scratch, COMMAND_RUNNER_FILENAME)
        runnerFile.text = task.processor.normalizeScript(runnerText)
        runnerFile.setExecutable(true)

        /*
         * save the reference to the scriptFile
         */
        task.script = scriptFile


        /*
         * Launch the job using the DRMAA api
         */
        def session = factory.getSession()
        session.init("")
        JobTemplate job = session.createJobTemplate()
        job.setWorkingDirectory( scratch.absolutePath )
        job.setRemoteCommand( runnerFile.absolutePath )
        job.setJobEnvironment( task.processor.getProcessEnvironment() )
        job.setJobName("nf-${task.processor.name}-${task.index}")
        job.setJoinFiles(true)
        job.setOutputPath('/dev/null')

        if( taskConfig['nativeGridOptions'] ) {
            def val = taskConfig['nativeGridOptions']
            if( val instanceof Collection ) {
                job.setNativeSpecification( val.join(' ') )
            }
            else {
                job.setNativeSpecification( val.toString() )
            }
        }

        def jobId = session.runJob(job);
        task.status = TaskRun.Status.RUNNING
        log.debug "DRMAA job submitted with id: $jobId"



        // -- save the 'sub' process output
        def subOutFile = new File(scratch, JOB_OUT_FILENAME)
        def subOutStream = new BufferedOutputStream(new FileOutputStream(subOutFile))
        ByteDumper subDumper = new ByteDumper(process.getInputStream(), {  byte[] data, int len -> subOutStream.write(data,0,len) } )
        subDumper.setName("sub_${task.name}")
        subDumper.start()

        try {
            // -- wait the the process completes
            task.exitCode = process.waitFor()
            def success = task.exitCode in taskConfig.validExitCodes
            log.debug "Task completeted > ${task.name} -- exit code: ${task.exitCode}; accepted code(s): ${taskConfig.validExitCodes.join(',')}"

            subDumper.await(500)
            subOutStream.close()

            // there may be very loooong delay over NFS, wait at least one minute
            if( success ) {
                Duration.waitFor('60s') { cmdOutFile.exists() }
            }
            if( cmdOutFile.exists() && taskConfig.echo ) {
                print cmdOutFile.text
            }

        }
        finally {
            subDumper.terminate()

            // make sure to release all resources
            IOUtils.closeQuietly(process.in)
            IOUtils.closeQuietly(process.out)
            IOUtils.closeQuietly(process.err)
            process.destroy()

            task.output = collectResultFile(task,'-')
        }


    }

    @Override
    getStdOutFile(TaskRun task) {
        return null  //To change body of implemented methods use File | Settings | File Templates.
    }
}
