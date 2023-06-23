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

package nextflow.processor

import java.nio.file.FileSystems
import java.nio.file.Files
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.Executor
import nextflow.file.FileHelper
import nextflow.util.CacheHelper
import nextflow.util.Escape

/**
 * Collect tasks and submit them as task groups to the underlying
 * executor.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskGroupCollector {

    private Executor executor

    private int groupSize

    private Lock sync = new ReentrantLock()

    private List<TaskHandler> group

    private boolean closed = false

    TaskGroupCollector(Executor executor, int groupSize) {
        this.executor = executor
        this.groupSize = groupSize
        this.group = new ArrayList<>(groupSize)
    }

    /**
     * Add a task to the current group, and submit the group when it
     * reaches the desired size.
     *
     * @param task
     */
    void collect(TaskRun task) {
        sync.lock()

        try {
            // submit task directly if the collector is closed
            // or if the task is retried (since it might have dynamic resources)
            if( closed || task.config.getAttempt() > 1 ) {
                executor.submit(task)
                return
            }

            // create task handler
            final handler = executor.createTaskHandler(task)

            // add task to the group
            group << handler

            // submit task group when it is ready
            if( group.size() == groupSize ) {
                submit0(group)
                group = new ArrayList<>(groupSize)
            }
        }
        finally {
            sync.unlock()
        }
    }

    /**
     * Close the collector, submitting any remaining tasks as a partial task group.
     */
    void close() {
        sync.lock()

        try {
            if( group.size() > 0 )
                submit0(group)

            closed = true
        }
        finally {
            sync.unlock()
        }
    }

    protected void submit0(List<TaskHandler> group) {
        // prepare child job launcher scripts
        for( TaskHandler handler : group )
            handler.prepareLauncher()

        // submit task group to the underlying executor
        executor.submit(createTaskGroup(group))
    }

    /**
     * Create the task run for a task group.
     *
     * @param group
     */
    protected TaskRun createTaskGroup(List<TaskHandler> group) {
        final tasks = group.collect( h -> h.task )
        final first = tasks.first()

        // compute hash and work directory
        final hash = CacheHelper.hasher( tasks.collect( t -> t.getHash().asLong() ) ).hash()
        final workDir = FileHelper.getWorkFolder(executor.getWorkDir(), hash)

        Files.createDirectories(workDir)

        // create wrapper script
        final script = createTaskGroupScript(group)

        // create task group
        return new TaskGroup(
            id: first.id,
            index: first.index,
            processor: first.processor,
            type: first.type,
            config: first.processor.config.createTaskConfig(),
            context: new TaskContext(first.processor),
            hash: hash,
            workDir: workDir,
            script: script,
            children: group
        )
    }

    /**
     * Create the wrapper script for a task group.
     *
     * @param group
     */
    protected String createTaskGroupScript(List<TaskHandler> group) {
        // get work directory and launch command for each task
        final workDirs = group.collect( h -> h.getWorkDir() )
        final args = group.first().getLaunchCommand().toArray() as String[]
        final cmd = Escape.cli(args).replaceAll(workDirs.first(), '\\${task_dir}')

        // create wrapper script
        """
        array=( ${workDirs.collect( p -> Escape.path(p) ).join(' ')} )
        for task_dir in \${array[@]}; do
            export task_dir
            ${cmd} || true
        done
        """.stripIndent().leftTrim()
    }

}
