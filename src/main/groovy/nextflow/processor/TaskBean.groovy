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

package nextflow.processor

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.executor.BashWrapperBuilder
import nextflow.util.MemoryUnit

/**
 * Serializable task value object. Holds configuration values required to
 * launch the task execution
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TaskBean implements Serializable, Cloneable {

    String taskId

    String name

    def input

    def scratch

    Map<String,String> environment

    String headerScript

    String wrapperScript

    String dockerImage

    List<String> moduleNames

    Path workDir

    Path targetDir

    String script

    List<String> shell

    Map dockerConfig

    String dockerCpuset

    MemoryUnit dockerMemory

    boolean statsEnabled

    String beforeScript

    String afterScript

    boolean executable

    Map<String,Path> inputFiles

    List<String> outputFiles

    String unstageStrategy

    Path sharedDir

    Path binDir

    @PackageScope
    TaskBean() {
        shell = BashWrapperBuilder.BASH
        inputFiles = [:]
        outputFiles = []
    }

    TaskBean(TaskRun task) {

        this.taskId = task.id
        this.name = task.name

        // set the input (when available)
        this.input = task.stdin
        this.scratch = task.scratch
        this.workDir = task.workDir
        this.targetDir = task.targetDir

        // set the environment
        // note: create a copy of the process environment to avoid concurrent
        // process executions override each others
        this.environment = new HashMap( task.getProcessor().getProcessEnvironment() )
        this.environment.putAll( task.getInputEnvironment() )

        this.moduleNames = task.config.getModule()
        this.shell = task.config.getShell() ?: BashWrapperBuilder.BASH
        this.script = task.getScript()
        this.beforeScript = task.config.beforeScript
        this.afterScript = task.config.afterScript

        // docker config
        this.dockerImage = task.getContainer()
        this.dockerConfig = task.getDockerConfig()
        this.dockerMemory = task.config.getMemory()
        this.executable = task.isContainerExecutable()

        // stats
        this.statsEnabled = task.getProcessor().getSession().statsEnabled

        this.inputFiles = task.getInputFilesMap()
        this.outputFiles = task.getOutputFilesNames()
        this.sharedDir = task.getProcessor().getSession().getWorkDir()
        this.binDir = task.getProcessor().getSession().getBinDir()
        this.unstageStrategy = task.getProcessor().getConfig().unstageStrategy

    }

    @Override
    TaskBean clone() {
        (TaskBean)super.clone()
    }

}
