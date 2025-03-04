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
import nextflow.processor.TaskRun
import nextflow.fusion.FusionHelper
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
        lines << ''
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
        // better handled as cluster options.
        // That being said, I'll preserve here the sorts of things I was thinking of for UWisc systems.
        // def requirements = []
        // def rank = []
        //     rank << "(TARGET.HasRotationalScratch == false)"
        // requirements << "(OpSys != WINDOWS)"

        result << "log = ${task.getWorkDirStr()}/${task.CMD_LOG}".toString()

        if ( ! isFusionEnabled() ) {
            // result << "out = ${task.getWorkDirStr()}/${task.CMD_OUTFILE}".toString()
            // result << "error = ${task.getWorkDirStr()}/${task.CMD_ERRFILE}".toString()
            result << "out = ${task.CMD_OUTFILE}".toString()
            result << "error = ${task.CMD_ERRFILE}".toString()
            result << "stream_out = true"
            result << "stream_error = true"
            result << "executable = ${task.CMD_RUN}".toString()
            if( task.isContainerEnabled() ) {
                result << "universe = container"
                result << "container_image = ${task.getContainer()}".toString()
            } else {
                result << "universe = vanilla"
            }

            // result << "transfer_files = NO" // note: this will result in jobs only being run in shared file systems, as God and Nextflow intended. HT Condor will only work with Nextflow and a shared filesystem, either a physical one (this case) or one provided by Fusion (in which case, Fusion's s3 filesystem will prov)
        } else {
            result << "out = ${task.getWorkDirStr()}/${task.CMD_OUTFILE}".toString()
            result << "error = ${task.getWorkDirStr()}/${task.CMD_ERRFILE}".toString()
            result << "stream_out = true"
            result << "stream_error = true"
            result << "getenv = true"
            result << "universe = vanilla"
            // executable will be added to manifest by CondorTaskHandler in Fusion setups
        }
        result << "transfer_executable = False" // handled by nextflow
        result << "transfer_output_files=\"\""  // ditto

            // result << "initialdir = ${task.getWorkDirStr()}".toString()

        if( task.config.getCpus()>1 ) {
            result << "request_cpus = ${task.config.getCpus()}".toString()
            result << "machine_count = 1"
        }

        if( task.config.getMemory() ) {
            result << "request_memory = ${task.config.getMemory()}".toString()
        }

        if( task.config.getDisk() ) {
            result << "request_disk = ${task.config.getDisk()}".toString()
        } else {
            result << "request_disk = 1 GB" // need minimum 1GB to run on CHTC servers
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

        result << "queue"
        result << ''
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
        ["condor_q", "-nobatch"]
    }


    static protected Map<String,QueueStatus> DECODE_STATUS = [
            'U': QueueStatus.PENDING,   // Unexpanded
            'I': QueueStatus.PENDING,   // Idle
            'R': QueueStatus.RUNNING,   // Running
            'X': QueueStatus.ERROR,     // Removed
            'C': QueueStatus.DONE,      // Completed
            'H': QueueStatus.HOLD,      // Held
            'E': QueueStatus.ERROR      // Error
    ]


    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text) {
        final result = new LinkedHashMap<String, QueueStatus>()
        if( !text ) return result

        boolean started = false
        def itr = text.readLines().iterator()
        while( itr.hasNext() ) {
            String line = itr.next()
            if( !started ) {
                started = line.startsWith(' ID ')
                continue
            }

            if( !line.trim() ) {
                break
            }

            def cols = line.tokenize(' ')
            def id = cols[0]
            def st = cols[5]
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




    @InheritConstructors
    static class CondorWrapperBuilder extends BashWrapperBuilder {

        String manifest // This contains the job description for condor. Unlike most systems, Condor will handle running the apptainer command; we specify the environment and command run in the manifest.

        Path build() {
            final wrapper = super.build()
            // give execute permission to wrapper file
            wrapper.setExecutable(true)
            // save the condor manifest
            this.workDir.resolve(CMD_CONDOR).text = manifest
            return wrapper
        }
    }

// ## What I think I want to do here:
// Create and submit a submit file that runs the fusion command on a container universe with args and environment specified as fusion wants

    @CompileStatic
    class CondorTaskHandler extends GridTaskHandler {
        // Prior code, and context this will be launched in:
        protected BashWrapperBuilder createTaskWrapper(TaskRun task) {
            return fusionEnabled()
                ? fusionLauncher()
                : executor.createBashWrapperBuilder(task)
        }

        protected String stdinLauncherScript() {
            return fusionEnabled() ? fusionStdinWrapper() : wrapperFile.text
        }

        // protected String fusionStdinWrapper() {
        //     // I think the issue boils down to a turf dispute.
        //     // Nextflow wants to operate in vanilla universe and execute the apptainer command.
        //     // Condor wants to be given the command to execute in the container, and generate that command on its machine.
        //     // Because Condor runs on server-controlled environments, I think it makes sense to yield command generation to condor.
        //     // I'm just not sure it is kosher to run "apptainer xyz" in a vanilla universe environment w/ Condor.
        //     // If it is, then we should let nextflow generate the apptainer command.
        //     // But I'm not sure it is ok to do that. I think that is the source of our problems.
        //     final submit = fusionSubmitCli()
        //     final launcher = fusionLauncher()
        //     final config = task.getContainerConfig()
        //     final containerOpts = task.config.getContainerOptions()
        //     final fusion_generated_command = FusionHelper.runWithContainer(launcher, config, task.getContainer(), containerOpts, submit)

        //     def condor_submit_commands = []
        //     // result << "container_image = ${task.getContainer()}".toString()
        //     //  containerConfig.getEnvWhitelist() +
        //     def environment = [:]
        //     environment += launcher.fusionEnv()
        //     environment = environment.each { key, val -> "${key}=${val}" }.collect().join(' ').toString()
        //     environment += task.getContainerConfig().getEnvWhitelist().join(' ')
        //     condor_submit_commands << "env = \"" + environment + "\""

        //     def submit_list = fusionSubmitCli()
        //     def executable = submit_list[0]
        //     def arguments = submit_list.drop(1).join(' ')

        //     condor_submit_commands << "executable = ${executable}".toString()
        //     condor_submit_commands << "arguments = ${arguments}".toString()
        //     // result << "container_target_dir = "
        //     // result << "container_options = ${options}" // simply do not get to specify these in htcondor
        //     condor_submit_commands = condor_submit_commands.join(' ')
        //     println(condor_submit_commands + '\n' + submitDirective(task))

        //     // result << "executable = /usr/bin/fusion"
        //     // result << "arguments =  bash \'${task.getWorkDirStr()}/.command.run\'"

        //     // create an inline script to launch the job execution
        //     return submitDirective(task) + cmd + '\n'
        // }

        protected String fusionStdinWrapper() {
            // goal: to transform this: 'set +u; env - PATH="$PATH" ${TMP:+APPTAINERENV_TMP="$TMP"} ${TMPDIR:+APPTAINERENV_TMPDIR="$TMPDIR"} APPTAINERENV_FUSION_WORK="/fusion/s3/dwerling-bucket-01/work/0d/c183147541efe18db4986e879aa423" APPTAINERENV_FUSION_LICENSE_TOKEN="XXX" APPTAINERENV_AWS_S3_ENDPOINT="https://campus.s3.wisc.edu:443" APPTAINERENV_AWS_ACCESS_KEY_ID="Xg1fge9p5Zp6nBB90fKS" APPTAINERENV_FUSION_TAGS="[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)" APPTAINERENV_AWS_SECRET_ACCESS_KEY="XXX" apptainer exec --no-home --pid /home/lalli/chtc_nf/.nextflow/singularity/wave.seqera.io-wt-095ce0d9cb72-nextflow-bash-latest.img /usr/bin/fusion bash '/fusion/s3/dwerling-bucket-01/work/0d/c183147541efe18db4986e879aa423/.command.run'
            // to condor submit file speak.
            // Condor submit file wants three things:
            // because the executable (apptainer) is already on the machine, I need to specify transfer_executable=False
            FusionHelper.runWithContainer(fusionLauncher(), task.getContainerConfig(), task.getContainer(), task.config.getContainerOptions(), fusionSubmitCli())
            final submit = fusionSubmitCli()
            final launcher = fusionLauncher()
            final config = task.getContainerConfig()
            final containerConfig = task.getContainerConfig()
            final containerOpts = task.config.getContainerOptions()
            final cmd = FusionHelper.runWithContainer(launcher, config, task.getContainer(), containerOpts, submit)
            // println ('launcher:\n')
            // println (fusionLauncher())
            // println ('task.getContainerConfig():\n')
            // println (task.getContainerConfig())
            // println ('task.getContainer():\n')
            // println (task.getContainer())
            // println ('task.config.getContainerOptions():\n')
            // println (containerOpts)
            // println ('fusionSubmitCli():\n')
            // println (submit)
            // println ('FusionHelper.runWithContainer(launcher, config, task.getContainer(), containerOpts, submit):')
            // println (cmd)
            // println ("containerConfig.getEnvWhitelist():")
            // println (containerConfig.getEnvWhitelist())
            // println ("task.getContainerConfig().getEnvWhitelist():")
            // println (task.getContainerConfig().getEnvWhitelist())
            // println (task.getContainerConfig().getEnvWhitelist().getClass())
            // println ("launcher.fusionEnv()")
            // println (launcher.fusionEnv())
            // println (launcher.fusionEnv()[1])
            // println (launcher.fusionEnv().getClass())
            // println ("containerConfig.runOptions as String")
            // println (containerConfig.runOptions as String)
            // println ("containerConfig.fusionOptions()")
            // println (containerConfig.fusionOptions())


            // create an inline script to launch the job execution
            // println('#!/bin/bash\n' + submitDirective(task) + cmd + '\n')
            // return List.of('condor_submit', '.condor.submit', '-terse', '-queue', '1')
            // println( '#!/bin/bash\n' + submitDirective(task) + cmd + '\n')
            // return '#!/bin/bash\n' + submitDirective(task) + cmd + '\n'

            def result = []
            result << "transfer_executable = false"
            // result << "container_image = ${task.getContainer()}".toString()
            //  containerConfig.getEnvWhitelist() +
            // def environment = [:]
            def environment = launcher.fusionEnv()
            environment = environment.each { key, val -> "${key}=${val}" }.collect()
            environment += task.getContainerConfig().getEnvWhitelist()

            result << "env = \"" + environment.join(' ').toString() + "\""

            def executable = submit[0]
            def arguments = submit.drop(1).join(' ')

            result << "executable = ${executable}".toString()
            result << "arguments = ${arguments}".toString()
            // result << "container_target_dir = "
            // result << "container_options = ${options}" // simply do not get to specify these in htcondor
            result = result.join('\n')
            println(submitDirective(task))
            println(result + '\nqueue\n')
            return result + '\nqueue\n'

            // create an inline script to launch the job execution
            // println('#!/bin/bash\n' + submitDirective(task) + cmd + '\n')
            // return List.of('condor_submit', '.condor.submit', '-terse', '-queue', '1')
            // println( '#!/bin/bash\n' + submitDirective(task) + cmd + '\n')
            // return '#!/bin/bash\n' + submitDirective(task) + cmd + '\n'
        }
    }
    // /*
    //  * Prepare and launch the task in the underlying execution platform
    //  */
    // @Override
    // GridTaskHandler createTaskHandler(TaskRun task) {
    //     assert task
    //     assert task.workDir

    //     new CondorTaskHandler(task)//, this)
    // }
}