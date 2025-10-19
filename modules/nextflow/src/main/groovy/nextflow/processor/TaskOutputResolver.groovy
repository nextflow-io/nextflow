/*
 * Copyright 2013-2025, Seqera Labs
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

import java.nio.file.Path

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.exception.IllegalArityException
import nextflow.exception.MissingFileException
import nextflow.exception.MissingValueException
import nextflow.script.params.v2.ProcessFileOutput
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * Implements the resolution of task outputs.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskOutputResolver implements Map<String,Object> {

    private Map<String,ProcessFileOutput> declaredFiles

    private TaskRun task

    @Delegate
    private Map<String,Object> delegate

    TaskOutputResolver(Map<String,ProcessFileOutput> declaredFiles, TaskRun task) {
        this.declaredFiles = declaredFiles
        this.task = task
        this.delegate = task.context
    }

    /**
     * Get an environment variable from the task environment.
     *
     * The underscore is needed to prevent calls from being dispatched
     * to Nextflow.env().
     *
     * @param name
     */
    String _env(String name) {
        final result = env0(null).get(name)

        if( result == null )
            throw new MissingValueException("Missing environment variable: $name")

        return result
    }

    /**
     * Get the result of an eval command from the task environment.
     *
     * @param name
     */
    String eval(String name) {
        final evalCmds = task.getOutputEvals()
        final result = env0(evalCmds).get(name)

        if( result == null )
            throw new MissingValueException("Missing result of eval command: '${evalCmds.get(name)}'")

        return result
    }

    @Memoized
    private Map env0(Map<String,String> evalCmds) {
        new TaskEnvCollector(task.workDir, evalCmds).collect()
    }

    /**
     * Get a file from the task environment.
     *
     * The underscore is needed to prevent calls from being dispatched
     * to Nextflow.file().
     *
     * @param key
     */
    Path _file(Map opts=[:], String key) {
        final param = declaredFiles.get(key)
        final filePattern = param.getFilePattern(delegate)
        if( filePattern.startsWith('/') )
            throw new IllegalArgumentException("Process output file '${filePattern}' in `${task.lazyName()}` is an absolute path")

        final allFiles = files0(filePattern, opts)
        if( allFiles.isEmpty() && opts.optional )
            return null
        if( allFiles.size() != 1 )
            throw new IllegalArityException("Process output file '${filePattern}' in `${task.lazyName()}` yielded ${allFiles.size()} files but expected only one")

        final result = allFiles.first()
        task.outputFiles.add(result)
        return result
    }

    /**
     * Get a collection of files from the task environment.
     *
     * The underscore is needed to prevent calls from being dispatched
     * to Nextflow.files().
     *
     * @param key
     */
    Set<Path> _files(Map opts=[:], String key) {
        final param = declaredFiles.get(key)
        final filePattern = param.getFilePattern(delegate)
        if( filePattern.startsWith('/') )
            throw new IllegalArgumentException("Process output glob '${filePattern}' in `${task.lazyName()}` is an absolute path")

        final allFiles = files0(filePattern, opts)
        task.outputFiles.addAll(allFiles)
        return allFiles.toSet()
    }

    @Memoized
    private List<Path> files0(String filePattern, Map opts) {
        new TaskFileCollector([filePattern], opts, task).collect()
    }

    /**
     * Get the standard output from the task environment.
     */
    Object stdout() {
        final value = task.@stdout

        if( value == null )
            throw new IllegalArgumentException("Missing 'stdout' for process > ${task.lazyName()}")

        if( value instanceof Path && !value.exists() )
            throw new MissingFileException("Missing 'stdout' file: ${value.toUriString()} for process > ${task.lazyName()}")

        return value instanceof Path ? value.text : value?.toString()
    }

    /**
     * Get a variable from the task context.
     *
     * @param name
     */
    @Override
    @CompileDynamic
    Object get(Object name) {
        try {
            return InvokerHelper.getProperty(delegate, name)
        }
        catch( MissingPropertyException e ) {
            throw new MissingValueException("Missing variable in process output: ${e.property}")
        }
    }
}
