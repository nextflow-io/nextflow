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
import nextflow.executor.TaskArrayExecutor
import nextflow.file.FileHelper
import nextflow.util.CacheHelper
import nextflow.util.Escape
/**
 * Task monitor that batches tasks and submits them as job arrays
 * to an underlying task monitor.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskArrayCollector {

    /**
     * The set of directives which are used by the job array.
     */
    private static final List<String> ARRAY_DIRECTIVES = [
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

    private TaskArrayExecutor executor

    private int arraySize

    private Lock sync = new ReentrantLock()

    private List<TaskRun> array

    private boolean closed = false

    TaskArrayCollector(TaskProcessor processor, Executor executor, int arraySize) {
        if( executor !instanceof TaskArrayExecutor )
            throw new IllegalArgumentException("Executor '${executor.name}' does not support job arrays")

        this.processor = processor
        this.executor = (TaskArrayExecutor)executor
        this.arraySize = arraySize
        this.array = new ArrayList<>(arraySize)
    }

    /**
     * Add a task to the current array, and submit the array when it
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

            // add task to the array
            array.add(task)

            // submit job array when it is ready
            if( array.size() == arraySize ) {
                executor.submit(createTaskArray(array))
                array = new ArrayList<>(arraySize)
            }
        }
        finally {
            sync.unlock()
        }
    }

    /**
     * Close the collector, submitting any remaining tasks as a partial job array.
     */
    void close() {
        sync.lock()
        try {
            if( array.size() == 1 ) {
                executor.submit(array.first())
            }
            else if( array.size() > 0 ) {
                executor.submit(createTaskArray(array))
                array = null
            }
            closed = true
        }
        finally {
            sync.unlock()
        }
    }

    /**
     * Create the task run for a job array.
     *
     * @param tasks
     */
    protected TaskArrayRun createTaskArray(List<TaskRun> tasks) {
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
        final script = createArrayTaskScript(handlers)
        log.debug "Creating task array run >> $workDir\n$script"
        
        // create config for job array
        final rawConfig = new HashMap<String,Object>(ARRAY_DIRECTIVES.size())
        for( final key : ARRAY_DIRECTIVES ) {
            final value = processor.config.get(key)
            if( value != null )
                rawConfig[key] = value
        }

        // create job array
        final first = tasks.min( t -> t.index )
        final taskArray = new TaskArrayRun(
            id: first.id,
            index: first.index,
            processor: processor,
            type: processor.taskBody.type,
            config: new TaskConfig(rawConfig),
            context: new TaskContext(processor),
            hash: hash,
            workDir: workDir,
            script: script,
            children: handlers
        )
        taskArray.config.context = taskArray.context
        taskArray.config.process = taskArray.processor.name
        taskArray.config.executor = taskArray.processor.executor.name

        return taskArray
    }

    /**
     * Create the wrapper script for a job array.
     *
     * @param array
     */
    protected String createArrayTaskScript(List<TaskHandler> array) {
        // get work directory and launch command for each task
        final workDirs = array.collect( h -> executor.getArrayWorkDir(h) )
        """
        array=( ${workDirs.collect( p -> Escape.path(p) ).join(' ')} )
        export nxf_array_task_dir=${getArrayIndexRef()}
        ${executor.getArrayLaunchCommand('$nxf_array_task_dir')}
        """.stripIndent().leftTrim()
    }

    protected String getArrayIndexRef() {
        final name = executor.getArrayIndexName()
        final start = executor.getArrayIndexStart()
        final index = start > 0 ? "${name} - ${start}" : name
        return '${array[' + index + ']}'
    }

}
