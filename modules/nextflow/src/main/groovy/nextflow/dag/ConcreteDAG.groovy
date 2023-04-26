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

package nextflow.dag

import java.nio.file.Path
import java.util.regex.Pattern
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.MapConstructor
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
import nextflow.script.params.FileOutParam
/**
 * Model the conrete (task) graph of a pipeline execution.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
class ConcreteDAG {

    private Lock sync = new ReentrantLock()

    private Map<String,Task> nodes = new HashMap<>(100)

    Map<String,Task> getNodes() {
        nodes
    }

    /**
     * Add a task to the graph
     *
     * @param task
     */
    void addTask(TaskRun task) {
        final hash = task.hash.toString()
        final label = "[${hash.substring(0,2)}/${hash.substring(2,8)}] ${task.name}"
        final inputs = task.getInputFilesMap()
            .collect { name, path ->
                new Input(name: name, path: path, predecessor: getPredecessorHash(path))
            }

        sync.lock()
        try {
            nodes[hash] = new Task(
                index: nodes.size(),
                label: label,
                inputs: inputs
            )
        }
        finally {
            sync.unlock()
        }
    }

    static public String getPredecessorHash(Path path) {
        final pattern = Pattern.compile('.*/([0-9a-f]{2}/[0-9a-f]{30})')
        final matcher = pattern.matcher(path.toString())

        matcher.find() ? matcher.group(1).replace('/', '') : null
    }

    /**
     * Add a task's outputs to the graph
     *
     * @param task
     */
    void addTaskOutputs(TaskRun task) {
        final hash = task.hash.toString()
        final outputs = task.getOutputsByType(FileOutParam)
            .values()
            .flatten()
            .collect { path ->
                new Output(name: path.name, path: path)
            }

        sync.lock()
        try {
            nodes[hash].outputs = outputs
        }
        finally {
            sync.unlock()
        }
    }

    @MapConstructor
    @ToString(includeNames = true, includes = 'label', includePackage=false)
    static protected class Task {
        int index
        String label
        List<Input> inputs
        List<Output> outputs

        String getSlug() { "t${index}" }
    }

    @MapConstructor
    static protected class Input {
        String name
        Path path
        String predecessor
    }

    @MapConstructor
    static protected class Output {
        String name
        Path path
    }

}
