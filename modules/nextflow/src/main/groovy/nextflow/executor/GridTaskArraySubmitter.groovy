/*
 * Copyright 2013-2023, Seqera Labs
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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessFailedException
import nextflow.exception.ProcessNonZeroExitStatusException
import nextflow.fusion.FusionHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.util.CmdLineHelper
/**
 * Submit tasks as an array job for a grid executor.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class GridTaskArraySubmitter extends TaskArraySubmitter implements SubmitJobAware {

    private AbstractGridExecutor executor

    GridTaskArraySubmitter(List<TaskHandler> array, AbstractGridExecutor executor) {
        super(array)
        this.executor = executor
    }

    @Override
    AbstractGridExecutor getExecutor() { executor }

    @Override
    TaskRun getTask() { array.first().getTask() }

    @Override
    protected void submit() {
        ProcessBuilder builder = null
        try {
            // -- create the array job script
            final launcherScript = getLauncherScript()

            // -- create the submit command
            builder = createProcessBuilder(true)

            // -- submit the array job with a retryable strategy
            final result = safeExecute( () -> launchProcess(builder, launcherScript) )
            final jobId = (String)executor.parseJobId(result)

            // -- set the job id and status of each task
            this.setJobId(jobId)
            this.setStatus(TaskStatus.SUBMITTED)

            log.debug "[${executor.name.toUpperCase()}] submitted array job > jobId: ${jobId}"
        }
        catch( Exception e ) {
            // update task exit status and message
            if( e instanceof ProcessNonZeroExitStatusException ) {
                for( TaskHandler handler : array ) {
                    handler.task.exitStatus = e.getExitStatus()
                    handler.task.stdout = e.getReason()
                    handler.task.script = e.getCommand()
                }
            }
            else {
                for( TaskHandler handler : array )
                    handler.task.script = builder ? CmdLineHelper.toLine(builder.command()) : null
            }
            this.setStatus(TaskStatus.COMPLETED)
            throw new ProcessFailedException("Error submitting array job for execution", e)
        }
    }

    protected String getLauncherScript() {
        return fusionEnabled()
            ? fusionLauncherScript()
            : classicLauncherScript()
    }

    protected String fusionLauncherScript() {
        final remoteLog = task.workDir.resolve(TaskRun.CMD_LOG).toString()
        final fusionWorkDir = FusionHelper.toContainerMount(task.workDir).toString()
        final arrayHeaders = executor.getArrayHeaders(array.size(), getTask())
        final arrayIndexName = executor.getArrayIndexName()
        final workDirs = array
            .collect { handler -> FusionHelper.toContainerMount(handler.task.workDir) }
            .join(' ')

        final cmd = FusionHelper.runWithContainer(
            fusionLauncher(),
            task.getContainerConfig(),
            task.getContainer(),
            fusionSubmitCli() )

        final builder = new StringBuilder()
            << '#!/bin/bash\n'
            << arrayHeaders.replace(remoteLog, '/dev/null')
            << "declare -a array=( ${workDirs} )\n"
            << cmd.replace(fusionWorkDir, "\${array[\$${arrayIndexName}]}")

        return builder.toString()
    }

    protected String classicLauncherScript() {
        final arrayHeaders = executor.getArrayHeaders(array.size(), getTask())
        final arrayIndexName = executor.getArrayIndexName()
        final workDirs = array
            .collect { handler -> handler.task.workDir }
            .join(' ')

        final builder = new StringBuilder()
            << '#!/bin/bash\n'
            << arrayHeaders
            << "declare -a array=( ${workDirs} )\n"
            << "bash \${array[\$${arrayIndexName}]}/${TaskRun.CMD_RUN}\n"

        return builder.toString()
    }
 
    protected void setJobId(String jobId) {
        array.eachWithIndex { handler, i ->
            ((GridTaskHandler)handler).setJobId(executor.getArrayTaskId(jobId, i))
        }
    }

    protected void setStatus(TaskStatus status) {
        for( TaskHandler handler : array )
            handler.setStatus(status)
    }

}
