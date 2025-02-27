/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
import nextflow.fusion.FusionHelper
import nextflow.processor.TaskRun
/**
 * HTCondor executor
 *
 * See https://research.cs.wisc.edu/htcondor/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CondorExecutor extends AbstractGridExecutor {

    static final public String CMD_CONDOR = '.command.condor'

    final protected BashWrapperBuilder createBashWrapperBuilder(TaskRun task) {
        // creates the wrapper script
        final builder = new BashWrapperBuilder(task)
        builder.manifest = getDirectivesText(task)
        return builder
    }

    protected String getDirectivesText(TaskRun task) {
        def lines = getDirectives(task)
        lines << ''
        lines.join('\n')
    }


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

        result << "log = ${TaskRun.CMD_LOG}"

        if ( ! isFusionEnabled() ) {
            result << "out = ${TaskRun.CMD_OUTFILE}"
            result << "error = ${TaskRun.CMD_ERRFILE}"
            result << "stream_out = true"
            result << "stream_error = true"
            result << "executable = ${TaskRun.CMD_RUN}"
            // result << "transfer_files = NO" // note: this will result in jobs only being run in shared file systems, as God and Nextflow intended. HT Condor will only work with Nextflow and a shared filesystem, either a physical one (this case) or one provided by Fusion (in which case, Fusion's s3 filesystem will prov)
        } else {
            result << "getenv = true"
            // executable will be added by CondorTaskHandler in Fusion setups
        }
        result << "transfer_executable = False" // handled by nextflow
        result << "transfer_output_files=\"\""  // ditto

        if( task.isContainerEnabled() ) {
            result << "universe = container"
            result << "container_image = ${task.getContainer()}"
        } else {
            result << "universe = vanilla"
            result << "initialdir = ${TaskRun.workDir}"
        }

        if( task.config.getCpus() > 1 ) {
            result << "request_cpus = ${task.config.getCpus()}"
            result << "machine_count = 1"
        }

        if( task.config.getMemory() ) {
            result << "request_memory = ${task.config.getMemory()}"
        }

        if( task.config.getDisk() ) {
            result << "request_disk = ${task.config.getDisk()}"
        } else {
            result << "request_disk = 1 GB" // need minimum 1GB to run on CHTC servers
        }

        if( task.config.getTime() ) {
            result << "periodic_remove = (RemoteWallClockTime - CumulativeSuspensionTime) > ${task.config.getTime().toSeconds()}"
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

        // if (requirements.size() > 0){
        //     result << "requirements = \"" + requirements.join(' && ') + "\""
        // }
        // if (rank.size() > 0){
        //     result << "rank = \"" + rank.join(' && ') + "\""
        // }

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


    static protected Map DECODE_STATUS = [
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
        def result = [:]
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

    /*
     * Prepare and launch the task in the underlying execution platform
     */
    @Override
    GridTaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir

        new CondorTaskHandler(task, this)
    }




    @InheritConstructors
    @CompileStatic
    class CondorTaskHandler extends GridTaskHandler {
        // Prior code, and context this will be launched in:
        // protected BashWrapperBuilder createTaskWrapper(TaskRun task) {
        //     return fusionEnabled()
        //         ? fusionLauncher()
        //         : executor.createBashWrapperBuilder(task)
        // }

        // protected String stdinLauncherScript() {
        //     return fusionEnabled() ? fusionStdinWrapper() : wrapperFile.text
        // }

        // protected String fusionStdinWrapper() {
        //     final submit = fusionSubmitCli()
        //     final launcher = fusionLauncher()
        //     final config = task.getContainerConfig()
        //     final containerOpts = task.config.getContainerOptions()
        //     final cmd = FusionHelper.runWithContainer(launcher, config, task.getContainer(), containerOpts, submit)
        //     // create an inline script to launch the job execution
        //     return '#!/bin/bash\n' + submitDirective(task) + cmd + '\n'
        // }

        protected String stdinLauncherScript() {
            FusionHelper.runWithContainer(fusionLauncher(), task.getContainerConfig(), task.getContainer(), task.config.getContainerOptions(), fusionSubmitCli())
            final submit = fusionSubmitCli()
            final launcher = fusionLauncher()
            final config = task.getContainerConfig()
            final containerConfig = task.getContainerConfig()
            final containerOpts = task.config.getContainerOptions()
            final cmd = FusionHelper.runWithContainer(launcher, config, task.getContainer(), containerOpts, submit)
            println ('launcher:\n')
            println (fusionLauncher())
            println ('task.getContainerConfig():\n')
            println (task.getContainerConfig())
            println ('task.getContainer():\n')
            println (task.getContainer())
            println ('task.config.getContainerOptions():\n')
            println (containerOpts)
            println ('fusionSubmitCli():\n')
            println (submit)
            println ('FusionHelper.runWithContainer(launcher, config, task.getContainer(), containerOpts, submit):')
            println (cmd)
            println ("containerConfig.getEnvWhitelist():")
            println (containerConfig.getEnvWhitelist())
            println ("task.getContainerConfig().getEnvWhitelist():")
            println (task.getContainerConfig().getEnvWhitelist())
            println (task.getContainerConfig().getEnvWhitelist().getClass())
            println ("launcher.fusionEnv()")
            println (launcher.fusionEnv())
            println (launcher.fusionEnv()[1])
            println (launcher.fusionEnv().getClass())
            println ("containerConfig.runOptions as String")
            println (containerConfig.runOptions as String)
            println ("containerConfig.fusionOptions()")
            println (containerConfig.fusionOptions())


            // create an inline script to launch the job execution
            // println('#!/bin/bash\n' + submitDirective(task) + cmd + '\n')
            // return List.of('condor_submit', '.condor.submit', '-terse', '-queue', '1')
            println( '#!/bin/bash\n' + submitDirective(task) + cmd + '\n')
            // return '#!/bin/bash\n' + submitDirective(task) + cmd + '\n'

            def result = []

            result << "container_image = ${task.getContainer()}"
            //  containerConfig.getEnvWhitelist() +
            def environment = [:]
            environment += launcher.fusionEnv()
            environment = environment.each { key, val -> "${key}=${val}" }.collect().join(' ').toString()
            environment += task.getContainerConfig().getEnvWhitelist().join(' ')
            result << "env = \"" + environment + "\""

            def submit_list = fusionSubmitCli()
            def executable = submit_list[0]
            def arguments = submit_list.drop(1).join(' ')

            result << "executable = ${executable}"
            result << "arguments = ${arguments}"
            // result << "container_target_dir = "
            // result << "container_options = ${options}" // simply do not get to specify these in htcondor
            result = result.join('\n')
            println(result + '\n' + submitDirective(task))
            return result + '\n' + submitDirective(task)

            // create an inline script to launch the job execution
            // println('#!/bin/bash\n' + submitDirective(task) + cmd + '\n')
            // return List.of('condor_submit', '.condor.submit', '-terse', '-queue', '1')
            // println( '#!/bin/bash\n' + submitDirective(task) + cmd + '\n')
            // return '#!/bin/bash\n' + submitDirective(task) + cmd + '\n'
        }
    }
}