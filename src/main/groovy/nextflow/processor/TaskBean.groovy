/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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
import nextflow.container.ContainerConfig
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

    String name

    def input

    def scratch

    Map<String,String> environment

    String headerScript

    String wrapperScript

    String containerImage

    List<String> moduleNames

    Path workDir

    Path targetDir

    String script

    List<String> shell

    ContainerConfig containerConfig

    String containerCpuset

    MemoryUnit containerMemory

    Path containerMount

    boolean statsEnabled

    String beforeScript

    String afterScript

    boolean containerExecutable

    boolean containerNative

    boolean containerEnabled

    String containerOptions

    Map<String,Path> inputFiles

    List<String> outputFiles

    String stageInMode

    String stageOutMode

    Path sharedDir

    Path binDir

    def cleanup

    @PackageScope
    TaskBean() {
        shell = BashWrapperBuilder.BASH
        inputFiles = [:]
        outputFiles = []
    }

    TaskBean(TaskRun task) {

        this.name = task.name

        // set the input (when available)
        this.input = task.stdin
        this.scratch = task.scratch
        this.workDir = task.workDir
        this.targetDir = task.targetDir

        // set the environment
        this.environment = task.getEnvironment()

        this.moduleNames = task.config.getModule()
        this.shell = task.config.getShell() ?: BashWrapperBuilder.BASH
        this.script = task.getScript()
        this.beforeScript = task.config.beforeScript
        this.afterScript = task.config.afterScript
        this.cleanup = task.config.cleanup

        // container config
        this.containerImage = task.getContainer()
        this.containerConfig = task.getContainerConfig()
        this.containerMemory = task.config.getMemory()
        this.containerExecutable = task.isContainerExecutable()
        this.containerNative = task.isContainerNative()
        this.containerEnabled = task.isContainerEnabled()
        this.containerOptions = task.config.containerOptions

        // stats
        this.statsEnabled = task.getProcessor().getSession().statsEnabled

        this.inputFiles = task.getInputFilesMap()
        this.outputFiles = task.getOutputFilesNames()
        this.sharedDir = task.getProcessor().getSession().getWorkDir()
        this.binDir = task.getProcessor().getSession().getBinDir()
        this.stageInMode = task.getProcessor().getConfig().stageInMode
        this.stageOutMode = task.getProcessor().getConfig().stageOutMode

    }

    @Override
    TaskBean clone() {
        (TaskBean)super.clone()
    }

}
