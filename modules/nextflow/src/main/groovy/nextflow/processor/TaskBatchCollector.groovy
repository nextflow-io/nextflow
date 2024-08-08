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
class TaskBatchCollector extends TaskCollector {

    private TaskProcessor processor

    private boolean parallel

    TaskBatchCollector(TaskProcessor processor, Executor executor, int batchSize, boolean parallel) {
        super(executor, batchSize)
        this.processor = processor
        this.parallel = parallel
    }

    /**
     * Create the task run for a task batch.
     *
     * @param tasks
     */
    @Override
    protected TaskRun createAggregateTask(List<TaskRun> tasks) {
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

        // create task batch
        final first = tasks.min( t -> t.index )
        final taskBatch = new TaskBatchRun(
            id: first.id,
            index: first.index,
            processor: processor,
            type: processor.taskBody.type,
            config: createAggregateConfig(processor),
            context: new TaskContext(processor),
            hash: hash,
            workDir: workDir,
            script: script,
            children: handlers
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
