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

    /**
     * The set of directives which are used by the task batch.
     */
    private static final List<String> SUBMIT_DIRECTIVES = [
            'accelerator',
            'arch',
            'clusterOptions',
            'cpus',
            'disk',
            'machineType',
            'memory',
            'queue',
            'resourceLabels',
            'resourceLimits',
            'time',
            // only needed for container-native executors and/or Fusion
            'container',
            'containerOptions',
    ]

    private TaskProcessor processor

    private Executor executor

    private int batchSize

    private boolean parallel

    private Lock sync = new ReentrantLock()

    private List<TaskRun> batch

    private boolean closed = false

    TaskBatchCollector(TaskProcessor processor, Executor executor, int batchSize, boolean parallel) {
        this.processor = processor
        this.executor = executor
        this.batchSize = batchSize
        this.parallel = parallel
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

            // add task to the batch
            batch.add(task)

            // submit task batch when it is ready
            if( batch.size() == batchSize ) {
                executor.submit(createTaskBatch(batch))
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
            if( batch.size() == 1 ) {
                executor.submit(batch.first())
            }
            else if( batch.size() > 0 ) {
                executor.submit(createTaskBatch(batch))
                batch = null
            }
            closed = true
        }
        finally {
            sync.unlock()
        }
    }

    /**
     * Create the task run for a task batch.
     *
     * @param tasks
     */
    protected TaskRun createTaskBatch(List<TaskRun> tasks) {
        // prepare child job launcher scripts
        final handlers = tasks.collect( t -> executor.createTaskHandler(t) )
        for( TaskHandler handler : handlers ) {
            handler.prepareLauncher()
        }

        // create work directory
        final hash = CacheHelper.hasher( tasks.collect( t -> t.getHash().asLong() ) ).hash()
        final workDir = FileHelper.getWorkFolder(executor.getWorkDir(), hash)
        Files.createDirectories(workDir)

        // create wrapper script
        final script = createBatchTaskScript(handlers)
        log.debug "Creating task batch run >> $workDir\n$script"

        // create config for task batch
        final rawConfig = new HashMap<String,Object>(SUBMIT_DIRECTIVES.size())
        for( final key : SUBMIT_DIRECTIVES ) {
            final value = processor.config.get(key)
            if( value != null )
                rawConfig[key] = value
        }

        // create task batch
        final first = tasks.min( t -> t.index )
        final taskBatch = new TaskBatchRun(
            id: first.id,
            index: first.index,
            processor: processor,
            type: processor.taskBody.type,
            config: new TaskConfig(rawConfig),
            context: new TaskContext(processor),
            hash: hash,
            workDir: workDir,
            script: script,
            children: tasks
        )
        taskBatch.config.context = taskBatch.context
        taskBatch.config.process = taskBatch.processor.name
        taskBatch.config.executor = taskBatch.processor.executor.name

        return taskBatch
    }

    /**
     * Create the wrapper script for a task batch.
     *
     * @param batch
     */
    protected String createBatchTaskScript(List<TaskHandler> batch) {
        final workDirs = batch.collect( h -> Escape.path(executor.getChildWorkDir(h)) )
        """
        array=( ${workDirs.join(' ')} )
        for nxf_batch_task_dir in \${array[@]}; do
            export nxf_batch_task_dir
            ${executor.getChildLaunchCommand('$nxf_batch_task_dir')} || true${parallel ? ' &' : ''}
        done
        ${parallel ? 'wait' : ''}
        """.stripIndent().leftTrim()
    }

}
