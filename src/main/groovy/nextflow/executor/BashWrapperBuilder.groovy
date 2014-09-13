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
import java.nio.file.Path

import groovy.util.logging.Slf4j
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.DockerBuilder
/**
 * Builder to create the BASH script which is used to
 * wrap and launch the user task
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class BashWrapperBuilder {

    final private TaskRun task

    def scratch

    def input

    Map<String,String> environment

    String stagingScript

    String unstagingScript

    String wrapperScript

    String dockerImage

    List<String> moduleNames

    Path workDir

    Path targetDir

    String script

    def shell

    Map dockerConfig

    Path dockerMount

    String name

    private runWithDocker

    BashWrapperBuilder( TaskRun task ) {
        this.task = task
        this.name = "nxf-" + task.hash?.toString()?.substring(0,8)

        // set the input (when available)
        this.input = task.stdin
        this.scratch = task.scratch
        this.workDir = task.workDir
        this.targetDir = task.targetDir

        // set the environment
        this.environment = task.processor.getProcessEnvironment()
        this.environment.putAll( task.getInputEnvironment() )

        this.moduleNames = task.processor.taskConfig.getModule()
        this.shell = task.processor.taskConfig.shell
        this.script = task.script.toString()

        // docker config
        this.dockerConfig = task.processor?.session?.config?.docker
        this.dockerImage = task.container ?: dockerConfig?.image
    }

    BashWrapperBuilder( Map params ) {
        log.trace "Wrapper params: $params"

        task = null
        this.name = params.name
        this.shell = params.shell ?: TaskConfig.DEFAULT_SHELL
        this.script = params.script?.toString()
        this.input = params.input
        this.scratch = params.scratch
        this.workDir = params.workDir
        this.targetDir = params.targetDir
        this.environment = params.environment

        // docker config
        this.dockerConfig = params.dockerConfig
        this.dockerImage = params.container ?: dockerConfig?.image
        this.dockerMount = params.dockerMount
    }

    /**
     * @return The bash script fragment to change to the 'scratch' directory if it has been specified in the task configuration
     */
    protected String changeToScratchDirectory() {

        if( scratch == null || scratch == false ) {
            return null
        }

        /*
         * when 'scratch' is defined as a bool value
         * try to use the 'TMP' variable, if does not exist fallback to a tmp folder
         */
        if( scratch == true ) {
            return 'NXF_SCRATCH=${TMPDIR:-`mktemp -d`} && cd $NXF_SCRATCH'
        }

        // convert to string for safety
        final scratchStr = scratch.toString()

        // when it is defined by a variable, just use it
        if( scratchStr.startsWith('$') ) {
            return "NXF_SCRATCH=\${${scratchStr.substring(1)}:-`mktemp -d`} && cd \$NXF_SCRATCH"
        }

        if( scratchStr.toLowerCase() in ['ramdisk','ram-disk']) {
            return 'NXF_SCRATCH=$(mktemp -d -p /dev/shm/) && cd $NXF_SCRATCH'
        }


        return "NXF_SCRATCH=\$(mktemp -d -p $scratch) && cd \$NXF_SCRATCH"

    }

    /**
     * Build up the BASH wrapper script file which will launch the user provided script
     * @return The {@code Path} of the created wrapper script
     */

    Path build() {
        assert workDir, "Missing 'workDir' property in BashWrapperBuilder object"
        assert script, "Missing 'script' property in BashWrapperBuilder object"

        final scriptFile = workDir.resolve(TaskRun.CMD_SCRIPT)
        final inputFile = workDir.resolve(TaskRun.CMD_INFILE)
        final environmentFile = workDir.resolve(TaskRun.CMD_ENV)
        final startedFile = workDir.resolve(TaskRun.CMD_START)
        final outputFile = workDir.resolve(TaskRun.CMD_OUTFILE)
        final exitedFile = workDir.resolve(TaskRun.CMD_EXIT)
        final wrapperFile = workDir.resolve(TaskRun.CMD_RUN)

        // set true when running with docker
        runWithDocker = dockerImage && dockerConfig?.enabled?.toString() == 'true'

        /*
         * the script file
         */
        final taskScript = scriptFile.text = TaskProcessor.normalizeScript(script, shell)

        /*
         * fetch the script interpreter
         */
        final interpreter = TaskProcessor.fetchInterpreter(taskScript)

        /*
         * save the input when required
         */
        if( input != null ) {
            inputFile.text = input
        }

        /*
         * save the environment
         */
        if( environment ) {
            // create the *bash* environment script
            environmentFile.text = TaskProcessor.bashEnvironmentScript(environment)
        }

        /*
         * create a script wrapper which do the following
         * 1 - move the TMP directory provided by the sge/oge grid engine
         * 2 - pipe the input stream
         * 3 - launch the user script
         * 4 - un-stage e.g. copy back the result files to the working folder
         */

        final ENDL = '\n'
        def wrapper = new StringBuilder()
        wrapper << '#!/bin/bash -Eeu' << ENDL
        wrapper << 'trap on_exit 1 2 3 15 ERR TERM USR1 USR2' << ENDL
        wrapper << 'function on_exit() { local exit_status=${1:-$?}; printf $exit_status > ' << exitedFile.toString()
        wrapper << '; exit $exit_status; }' << ENDL
        wrapper << 'touch ' << startedFile.toString() << ENDL

        // source the environment
        if( !runWithDocker ) {
            wrapper << '[ -f '<< environmentFile.toString() << ' ]' << ' && source ' << environmentFile.toString() << ENDL
        }

        // when a module is defined, load it
        moduleNames?.each { String name ->
            wrapper << 'module load ' << name << ENDL
        }

        // whenever it has to change to the scratch directory
        def changeDir = changeToScratchDirectory()
        if( changeDir ) {
            wrapper << changeDir << ENDL
        }

        // staging input files when required
        if( stagingScript ) {
            wrapper << stagingScript << ENDL
        }

        // execute the command script
        wrapper << '( '

        // execute by invoking the command through a Docker container
        DockerBuilder docker = null
        if( runWithDocker ) {
            docker = new DockerBuilder(dockerImage)
            if( task ) {
                docker.addMountForInputs( task.getInputFiles() )
                      .addMount( task.processor.session.workDir )
                      .addMount( task.processor.session.binDir )
            }

            // set the name
            docker.setName(this.name)

            if( dockerMount )
                docker.addMount(dockerMount)

            // set the environment
            if( !environmentFile.empty() )
                docker.addEnv( environmentFile )

            // turn on container remove by default
            if( !dockerConfig.containsKey('remove') )
                dockerConfig.remove = true

            // set up run docker params
            docker.params(dockerConfig)

            // extra rule for the 'auto' temp dir temp dir
            def temp = dockerConfig.temp?.toString()
            if( temp == 'auto' || temp == 'true' ) {
                docker.setTemp( changeDir ? '$NXF_SCRATCH' : '$(mktemp -d)' )
            }

            wrapper << docker.build() << ' '
        }

        wrapper << interpreter << ' ' << scriptFile.toString()
        if( input != null ) wrapper << ' < ' << inputFile.toString()
        wrapper << ' ) &> ' << outputFile.toAbsolutePath() << ENDL

        if( docker?.removeCommand ) {
            // remove the container in this way because 'docker run --rm'  fail in some cases -- see https://groups.google.com/d/msg/docker-user/0Ayim0wv2Ls/-mZ-ymGwg8EJ
            wrapper << docker.removeCommand << ' &>/dev/null || true &' << ENDL
        }

        if( (changeDir || workDir != targetDir) && unstagingScript  ) {
            wrapper << unstagingScript << ENDL
        }

        wrapper << 'on_exit' << ENDL

        wrapperFile.text = wrapperScript = wrapper.toString()
        return wrapperFile
    }




}
