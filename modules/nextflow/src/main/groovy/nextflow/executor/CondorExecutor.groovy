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
import groovy.transform.InheritConstructors
import nextflow.container.ContainerBuilder
import nextflow.processor.TaskRun
import nextflow.fusion.FusionHelper
import nextflow.file.FileHelper
import nextflow.extension.FilesEx
import nextflow.exception.ProcessException
import nextflow.util.MemoryUnit
import nextflow.SysEnv

import static java.nio.file.StandardOpenOption.*

import java.nio.file.FileSystemException
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path



/**
 * HTCondor executor
 *
 * See https://research.cs.wisc.edu/htcondor/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CondorExecutor extends AbstractGridExecutor {

    static final public String CMD_CONDOR = '.command.condor'

    final protected BashWrapperBuilder createBashWrapperBuilder(TaskRun task) {
        // creates the wrapper script
        final builder = new CondorWrapperBuilder(task)
        builder.manifest = getDirectivesText(task)
        return builder
    }

    protected String getDirectivesText(TaskRun task) {
        def lines = getDirectives(task)
        lines.join('\n')
    }

    // @Override
    // protected String getHeaderToken() {
    //     throw new UnsupportedOperationException()
    // }


    // Condor does not require a special token or header
    protected String getHeaderToken() { return '' }

    /**
     * Defines the jobs directive headers
     *
     * @param task
     * @return A multi-line string containing the job directives
     */
    String getHeaders( TaskRun task ) {
        return getDirectivesText(task)
    }


	// We currently assume that the system has been configured so that
	// anyone (user) who can run an HTCondor job can also run docker.  It's
	// also apparently a security worry to run Docker as root, so let's not.
    // https://github.com/htcondor/htcondor/blob/02c6bc70543951cf9d352e3fbf3343a925a47e3a/src/condor_utils/docker-api.cpp#L99C1-L102C1

    @Override
    protected List<String> getDirectives(TaskRun task, List<String> result) {

        result << "universe = vanilla"

        result << "out = ${TaskRun.CMD_OUTFILE}".toString()
        result << "error = ${TaskRun.CMD_ERRFILE}".toString()
        result << "log = .condor_runlog.uuid-${session.uniqueId}.log".toString()
        result << "getenv = true"

        result << "transfer_executable = False" // handled by nextflow
        result << "transfer_output_files=\"\""  // ditto

        if( task.config.getCpus()>1 ) {
            result << "request_cpus = ${task.config.getCpus()}".toString()
            result << "machine_count = 1"
        }

        if( task.config.getMemory() ) {
            result << "request_memory = ${task.config.getMemory()}".toString()
        }

        if( task.config.getDisk() ) {
            result << "request_disk = ${task.config.getDisk()}".toString()
        }

        if( task.config.getTime() ) {
            result << "periodic_remove = (RemoteWallClockTime - CumulativeSuspensionTime) > ${task.config.getTime().toSeconds()}".toString()
        }

        if( task.config.getClusterOptions() ) {
            def opts = task.config.getClusterOptions()
            if( opts instanceof Collection ) {
                result.addAll(opts as Collection)
            }
            else {
                result.addAll( opts.toString().tokenize(';\n').collect{ it.trim() })
            }
        }
        if ( ! pipeLauncherScript() ) {
            if ( task.isContainerEnabled() ){

            } else {
                result << "executable = ${task.CMD_RUN}".toString()
                result << "environment = ${task.getEnvironment()}".toString()
            }
            // if not containerized, the executable is .command.run, which we will specify in the get directives portion.
            // result << "executable = placeholder"
            // result << "arguments = placeholder"
            // result << "environment = placeholder"
            result << 'queue'
        }
        return result
    }


    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile) {
        return pipeLauncherScript()
                ? List.of('condor_submit', '-', '-terse',)
                : List.of('condor_submit', '-terse', CMD_CONDOR)
    }

    @Override
    def parseJobId(String text) {
        text.tokenize(' -')[0]
    }

    @Override
    protected List<String> getKillCommand() {
        ['condor_rm']
    }

    @Override
    protected List<String> queueStatusCommand(Object queue) {
        ["condor_history", "-userlog", ".condor_runlog.uuid-${session.uniqueId}.log".toString(), "-wide","-af:j", "JobStatus"]
    }


    static protected Map<String,QueueStatus> DECODE_STATUS = [
            'U': QueueStatus.PENDING,   // Unexpanded
            'I': QueueStatus.PENDING,   // Idle
            'R': QueueStatus.RUNNING,   // Running
            'X': QueueStatus.ERROR,     // Removed
            'C': QueueStatus.DONE,      // Completed
            'H': QueueStatus.HOLD,      // Held
            'E': QueueStatus.ERROR,      // Error
            // numeric options
            '0': QueueStatus.PENDING,   // Unexpanded
            '1': QueueStatus.PENDING,   // Idle
            '2': QueueStatus.RUNNING,   // Running
            '3': QueueStatus.ERROR,     // Removed
            '4': QueueStatus.DONE,      // Completed
            '5': QueueStatus.HOLD,      // Held
            '6': QueueStatus.ERROR      // Error
    ]


    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {
        final result = new LinkedHashMap<String, QueueStatus>()
        if( !text ) {
            println("escaping because !text")
            return result
        }

        def itr = text.readLines().iterator()
        while( itr.hasNext() ) {
            String line = itr.next()
            if( line.startsWith(' ID ') ) continue

            if( !line.trim() ) {
                break
            }

            def cols = line.tokenize(' ')
            def id = cols[0]
            def st = cols[1]
            result[id] = DECODE_STATUS[st]
        }

        return result
    }

    @Override
    protected boolean pipeLauncherScript() {
        return isFusionEnabled()
    }

    @Override
    boolean isFusionEnabled() {
        return FusionHelper.isFusionEnabled(session)
    }


    /*
     * Prepare and launch the task in the underlying execution platform
     */
    CondorTaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir

        new CondorTaskHandler(task, this)
    }


    /**
    * Handles a job execution in the underlying grid platform
    */
    @CompileStatic
    @InheritConstructors
    class CondorTaskHandler extends GridTaskHandler {
        // creates text of bash wrapper executable file to run on remote server
        protected String generateFusionBashWrapperCommand() {
            final submit = fusionSubmitCli()
            final launcher = fusionLauncher()
            final config = task.getContainerConfig()
            final containerOpts = task.config.getContainerOptions()
            final cmd = FusionHelper.runWithContainer(launcher, config, task.getContainer(), containerOpts, submit)

            return '#!/bin/bash\n' + cmd + '\n'
        }

        // creates condor submit file that is fed to stdin. The submit file specifies an executable
        protected String fusionStdinWrapper() {
            final fusionBashWrapperText = generateFusionBashWrapperCommand()

            final String tmp_launch_script = ".condor.${task.id}.${task.hash}.sh"
            final Path executable_file_name = FileHelper.getLocalTempPath().resolve(tmp_launch_script)
            // save the bash command to a script on disk
            executable_file_name.text = fusionBashWrapperText

            FilesEx.setExecutable(executable_file_name, true)
            def submit_file_commands = getDirectives(task)
            submit_file_commands << "executable = ${executable_file_name}".toString()
            submit_file_commands << "transfer_executable = True"
            submit_file_commands << "queue"
            submit_file_commands << ""

            return submit_file_commands.join('\n')

        }

        private static MemoryUnit DEFAULT_STAGE_FILE_THRESHOLD = MemoryUnit.of('1 MB')
        private static int DEFAULT_WRITE_BACK_OFF_BASE = 3
        private static int DEFAULT_WRITE_BACK_OFF_DELAY = 250
        private static int DEFAULT_WRITE_MAX_ATTEMPTS = 5

        private MemoryUnit stageFileThreshold = SysEnv.get('NXF_WRAPPER_STAGE_FILE_THRESHOLD') as MemoryUnit ?: DEFAULT_STAGE_FILE_THRESHOLD
        private int writeBackOffBase = SysEnv.get('NXF_WRAPPER_BACK_OFF_BASE') as Integer ?: DEFAULT_WRITE_BACK_OFF_BASE
        private int writeBackOffDelay = SysEnv.get('NXF_WRAPPER_BACK_OFF_DELAY') as Integer ?: DEFAULT_WRITE_BACK_OFF_DELAY
        private int writeMaxAttempts = SysEnv.get('NXF_WRAPPER_MAX_ATTEMPTS') as Integer ?: DEFAULT_WRITE_MAX_ATTEMPTS
        static protected boolean isRetryable0(Exception e) {
            if( e instanceof FileSystemException )
                return true
            if( e instanceof SocketException )
                return true
            if( e instanceof RuntimeException )
                return true
            if( e.class.getSimpleName() == 'HttpResponseException' )
                return true
            return false
    }

        private Path write0(Path path, String data) {
            int attempt=0
            while( true ) {
                try {
                    try (BufferedWriter writer=Files.newBufferedWriter(path, CREATE,WRITE,TRUNCATE_EXISTING)) {
                        writer.write(data)
                    }
                    return path
                }
                catch (Exception e) {
                    if( !this.isRetryable0(e) )
                        throw e
                    final isLocalFS = path.getFileSystem()==FileSystems.default
                    // the retry logic is needed for non-local file system such as S3.
                    // when the file is local fail without retrying
                    if( isLocalFS || ++attempt>=writeMaxAttempts )
                        throw new ProcessException("Unable to create file ${path.toUriString()}", e)
                    // use an exponential delay before making another attempt
                    final delay = (Math.pow(writeBackOffBase, attempt) as long) * writeBackOffDelay
                    Thread.sleep(delay)
                }
            }
        }
    }

    @InheritConstructors
    static class CondorWrapperBuilder extends BashWrapperBuilder {

        String manifest

        Path build() {
            final wrapper = super.build()
            // give execute permission to wrapper file
            wrapper.setExecutable(true)
            // save the condor manifest
            this.workDir.resolve(CMD_CONDOR).text = manifest
            return wrapper
        }


    }
}
