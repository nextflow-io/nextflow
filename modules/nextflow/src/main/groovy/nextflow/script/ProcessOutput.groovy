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

package nextflow.script

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.exception.MissingFileException
import nextflow.exception.MissingValueException
import nextflow.processor.TaskEnvCollector
import nextflow.processor.TaskFileCollecter
import nextflow.processor.TaskRun
/**
 * Models a process output.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ProcessOutput implements Cloneable {

    static enum Shortcuts { STDOUT }

    private ProcessOutputs parent

    private Object target

    private String name

    private boolean optional

    private DataflowWriteChannel channel

    ProcessOutput(ProcessOutputs parent, Object target, Map<String,?> opts) {
        this.parent = parent
        this.target = target

        if( opts.name )
            this.name = opts.name
        if( opts.optional )
            this.optional = true
    }

    String getName() {
        return name
    }

    void setChannel(DataflowWriteChannel channel) {
        this.channel = channel
    }

    DataflowWriteChannel getChannel() {
        return channel
    }

    Object resolve(TaskRun task) {
        final ctx = new ResolverContext(parent, optional, task)
        return resolve0(ctx)
    }

    private Object resolve0(ResolverContext ctx) {
        if( target == Shortcuts.STDOUT )
            return ctx.stdout()

        if( target instanceof Closure )
            return ctx.with(target)

        return target
    }

    static private class ResolverContext {

        private ProcessOutputs parent

        private boolean optional

        private TaskRun task

        ResolverContext(ProcessOutputs parent, boolean optional, TaskRun task) {
            this.parent = parent
            this.optional = optional
            this.task = task
        }

        /**
         * Get an environment variable from the task environment.
         *
         * @param key
         */
        String env(String key) {
            final varName = parent.env.get(key)
            final result = env0(task.workDir).get(varName)

            if( result == null && !optional )
                throw new MissingValueException("Missing environment variable: $varName")

            return result
        }

        @Memoized
        static private Map env0(Path workDir) {
            new TaskEnvCollector(workDir).collect()
        }

        /**
         * Get a file or list of files from the task environment.
         *
         * @param key
         */
        Object path(String key) {
            final param = parent.files.get(key)
            final result = new TaskFileCollecter(param, task).collect()

            if( result instanceof Path )
                task.outputFiles.add(result)
            else if( result instanceof Collection<Path> )
                task.outputFiles.addAll(result)

            return result
        }

        /**
         * Get the standard output from the task environment.
         */
        Object stdout() {
            final result = task.@stdout

            if( result == null && task.type == ScriptType.SCRIPTLET )
                throw new IllegalArgumentException("Missing 'stdout' for process > ${task.lazyName()}")

            if( result instanceof Path && !result.exists() )
                throw new MissingFileException("Missing 'stdout' file: ${result.toUriString()} for process > ${task.lazyName()}")

            return result
        }

        /**
         * Get a variable from the task context.
         *
         * @param name
         */
        @Override
        Object getProperty(String name) {
            if( name == 'stdout' )
                return stdout()

            try {
                return task.context.get(name)
            }
            catch( MissingPropertyException e ) {
                throw new MissingValueException("Missing variable in emit statement: ${e.property}")
            }
        }
    }
}
