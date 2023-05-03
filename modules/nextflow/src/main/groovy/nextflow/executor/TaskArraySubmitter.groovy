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

import java.nio.file.FileSystems
import java.nio.file.Files
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.FileHelper
import nextflow.fusion.FusionHelper
import nextflow.processor.TaskArray
import nextflow.processor.TaskContext
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.util.CacheHelper
import nextflow.util.Escape

/**
 * Submit tasks as an array job.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskArraySubmitter {

    private List<TaskHandler> array

    private TaskArrayAware executor

    private AtomicInteger collected = new AtomicInteger()

    TaskArraySubmitter(List<TaskHandler> array, TaskArrayAware executor) {
        this.array = array
        this.executor = executor
    }

    /**
     * Mark a task as ready to be submitted.
     *
     * When all tasks in the array are ready, the array job
     * will be submitted.
     *
     * @param handler
     */
    void collect(TaskHandler handler) {
        if( collected.incrementAndGet() == array.size() )
            submit()
    }

    /**
     * Submit the array job.
     */
    protected void submit() {
        final tasks = array.collect( h -> h.task )
        final first = tasks.first()

        // create work directory
        final hash = CacheHelper.hasher( tasks.collect( t -> t.getHash().asLong() ) ).hash()
        final workDir = FileHelper.getWorkFolder(executor.getWorkDir(), hash)

        Files.createDirectories(workDir)

        // create wrapper script
        final script = createWrapperScript(tasks)

        // create task handler
        final arrayTask = new TaskArray(
            id: first.id,
            index: first.index,
            processor: first.processor,
            type: first.type,
            config: first.processor.config.createTaskConfig(),
            context: new TaskContext(first.processor),
            hash: hash,
            workDir: workDir,
            script: script,
            children: array
        )
        final arrayHandler = executor.createTaskHandler(arrayTask)

        // submit array job
        arrayHandler.prepareLauncher()
        arrayHandler.submit()

        log.trace "Submitted array job ${arrayTask.name} > workDir: ${arrayTask.workDir}"
    }

    /**
     * Create the wrapper script for an array job.
     */
    protected String createWrapperScript(List<TaskRun> tasks) {
        // get work directory and launch command for each task
        def workDirs
        def cmd

        if( executor.workDir.fileSystem == FileSystems.default ) {
            workDirs = tasks.collect( t -> t.workDir.toString() )
            cmd = "cd \${task_dir} ; bash ${TaskRun.CMD_RUN} &> ${TaskRun.CMD_LOG}"
        }
        else {
            workDirs = executor.isFusionEnabled()
                ? tasks.collect( t -> FusionHelper.toContainerMount(t.workDir).toString() )
                : tasks.collect( t -> t.workDir.toUriString() )
            cmd = Escape.cli(array.first().getSubmitCommand().toArray() as String[])
            cmd = cmd.replaceAll(workDirs.first(), '\\${task_dir}')
        }

        // create wrapper script
        final arrayIndexName = executor.getArrayIndexName()

        """
        declare -a array=( ${workDirs.collect( p -> Escape.path(p) ).join(' ')} )
        export task_dir=\${array[\$${arrayIndexName}]}
        ${cmd}
        """.stripIndent().trim()
    }

}
