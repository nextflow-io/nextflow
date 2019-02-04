/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

    String containerImage

    Path condaEnv

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

    @Deprecated
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

        this.condaEnv = task.getCondaEnv()
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
