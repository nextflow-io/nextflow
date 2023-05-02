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
import nextflow.fusion.FusionHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
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
    String stdinLauncherScript() {
        def arrayDirective = executor.getArrayDirective(array.size(), task)
        def directives = executor.getDirectives(task, arrayDirective)
        def headers = executor.getHeaders(directives)
        def arrayIndexName = executor.getArrayIndexName()

        def workDirs
        def cmd

        if( fusionEnabled() ) {
            final remoteLog = task.workDir.resolve(TaskRun.CMD_LOG).toString()
            headers = headers.replaceAll(remoteLog, '/dev/null')

            workDirs = array.collect { h -> FusionHelper.toContainerMount(h.task.workDir).toString() }
            cmd = FusionHelper
                .runWithContainer(fusionLauncher(), task.getContainerConfig(), task.getContainer(), fusionSubmitCli())
                .replaceAll(workDirs.first(), '\\$task_dir')
        }
        else {
            workDirs = array.collect { h -> h.task.workDir.toString() }
            cmd = "bash \$task_dir/${TaskRun.CMD_RUN}"
        }

        final builder = new StringBuilder()
            << '#!/bin/bash\n'
            << headers
            << "declare -a array=( ${workDirs.join(' ')} )\n"
            << "task_dir=\${array[\$${arrayIndexName}]}\n"
            << "${cmd}\n"

        return builder.toString()
    }

    @Override
    void submit() {
        final jobId = submitJob(true)

        array.eachWithIndex { handler, i ->
            ((GridTaskHandler)handler).setJobId(executor.getArrayTaskId(jobId, i))
            handler.setStatus(TaskStatus.SUBMITTED)
        }

        log.debug "[${executor.name.toUpperCase()}] submitted array job > jobId: ${jobId}"
    }

    @Override
    Exception submitError(Exception e, String submitCommand) {
        for( TaskHandler handler : array )
            ((GridTaskHandler)handler).submitError(e, submitCommand)

        throw new ProcessFailedException("Error submitting array job for execution", e)
    }

}
