/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.script

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.TupleConstructor
import nextflow.Const
import nextflow.exception.ScriptRuntimeException
import nextflow.processor.PublishDir
import nextflow.processor.TaskRun
import nextflow.script.params.FileOutParam
import nextflow.script.ProcessConfig
/**
 * Models the workflow outputs definition and publishing
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class WorkflowPublisher {
    private List<PublisherEntry> publishers = []

    WorkflowPublisher(Path path, List<OutputCollection> collections) {
        for( def collection : collections ) {
            for( def selector : collection.selectors ) {
                final params = [
                    path: path.resolve(collection.path).resolve(selector.path),
                    pattern: selector.pattern,
                    failOnError: true
                ]
                publishers << new PublisherEntry(selector.name, PublishDir.create(params))
            }
        }
    }

    void publish(TaskRun task) {
        // collect task output files
        HashSet<Path> files = []
        final outputs = task.getOutputsByType(FileOutParam)
        for( Map.Entry entry : outputs ) {
            final value = entry.value
            if( value instanceof Path )
                files.add((Path)value)
            else if( value instanceof Collection<Path> )
                files.addAll(value)
            else if( value != null )
                throw new IllegalArgumentException("Unknown output file object [${value.class.name}]: ${value}")
        }

        // apply each publisher with matching process selector to task
        final processName = task.processor.name
        final simpleName = processName.split(Const.SCOPE_SEP).last()
        for( final entry : publishers ) {
            final selector = entry.selector
            final publisher = entry.publisher
            if( ProcessConfig.matchesSelector(simpleName, selector) || ProcessConfig.matchesSelector(processName, selector) )
                publisher.apply(files, task)
        }
    }

    @TupleConstructor
    private static class PublisherEntry {
        String selector
        PublishDir publisher
    }
}

@CompileStatic
@TupleConstructor
class OutputCollection {

    String path
    List<Selector> selectors
    Index index

    static class Selector {
        String name
        String path
        String pattern

        Selector(String name, String path, String pattern) {
            this.name = name
            this.path = path
            this.pattern = pattern
        }

        Selector(String name) {
            this(name, '.', null)
        }
    }

    @TupleConstructor
    static class Index {
        String format
        String path
    }

}
