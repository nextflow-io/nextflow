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
 * Collect tasks and submit them as task batches to the underlying
 * executor.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskBatchCollector {

    private Executor executor

    private int batchSize

    private Lock sync = new ReentrantLock()

    private List<TaskHandler> batch

    private boolean closed = false

    TaskBatchCollector(Executor executor, int batchSize) {
        this.executor = executor
        this.batchSize = batchSize
        this.batch = new ArrayList<>(batchSize)
    }

    /**
     * Add a task to the current batch, and submit the batch when it
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

            // add task to the batch
            batch << handler

            // submit task batch when it is ready
            if( batch.size() == batchSize ) {
                submit0(batch)
                batch = new ArrayList<>(batchSize)
            }
        }
        finally {
            sync.unlock()
        }
    }

    /**
     * Close the collector, submitting any remaining tasks as a partial task batch.
     */
    void close() {
        sync.lock()

        try {
            if( batch.size() > 0 )
                submit0(batch)

            closed = true
        }
        finally {
            sync.unlock()
        }
    }

    protected void submit0(List<TaskHandler> batch) {
        // prepare child job launcher scripts
        for( TaskHandler handler : batch )
            handler.prepareLauncher()

        // submit task batch to the underlying executor
        executor.submit(createTaskBatch(batch))
    }

    /**
     * Create the task run for a task batch.
     *
     * @param batch
     */
    protected TaskRun createTaskBatch(List<TaskHandler> batch) {
        final tasks = batch.collect( h -> h.task )
        final first = tasks.first()

        // compute hash and work directory
        final hash = CacheHelper.hasher( tasks.collect( t -> t.getHash().asLong() ) ).hash()
        final workDir = FileHelper.getWorkFolder(executor.getWorkDir(), hash)

        Files.createDirectories(workDir)

        // create wrapper script
        final script = createTaskBatchScript(batch)

        // create task batch
        return new TaskBatch(
            id: first.id,
            index: first.index,
            processor: first.processor,
            type: first.type,
            config: first.processor.config.createTaskConfig(),
            context: new TaskContext(first.processor),
            hash: hash,
            workDir: workDir,
            script: script,
            children: batch
        )
    }

    /**
     * Create the wrapper script for a task batch.
     *
     * @param batch
     */
    protected String createTaskBatchScript(List<TaskHandler> batch) {
        // get work directory and launch command for each task
        final workDirs = batch.collect( h -> h.getWorkDir() )
        final args = batch.first().getLaunchCommand().toArray() as String[]
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
