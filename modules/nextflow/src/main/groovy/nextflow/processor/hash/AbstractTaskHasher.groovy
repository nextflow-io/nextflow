/*
 * Copyright 2013-2026, Seqera Labs
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
package nextflow.processor.hash

import java.nio.file.Path

import com.google.common.hash.HashCode
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.UnexpectedException
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.CacheHelper
import nextflow.util.HashBuilder
/**
 * Common logic for task hash computation strategies.
 *
 * Subclasses implement {@link #createHashBuilder()} to configure
 * version-specific hashing behaviour.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class AbstractTaskHasher implements TaskHasher {

    protected TaskRun task

    protected TaskProcessor processor

    protected Session session

    AbstractTaskHasher(TaskRun task) {
        this.task = task
        this.processor = task.processor
        this.session = task.processor.session
    }

    /**
     * Create a {@link HashBuilder} configured for this hashing strategy.
     */
    abstract protected HashBuilder createHashBuilder()

    @Override
    HashCode compute() {

        final keys = new ArrayList<Object>()

        // add session UUID
        keys << session.uniqueId

        // add fully-qualified process name
        keys << task.processor.name

        // add source code of `script:` or `exec:` block
        //
        // - this allows task script to reference directives like `task.cpus`
        //   without invalidating the cache
        //
        // - references to local variables, global variables, and `task.ext`
        //   are included separately
        keys << task.source

        // add container fingerprint if present
        if( task.isContainerEnabled() )
            keys << task.getContainerFingerprint()

        // add the name and value of each task input
        for( final entry : task.inputs ) {
            keys.add( entry.key.name )
            keys.add( entry.value )
        }

        // add eval output commands
        final outEvals = task.getOutputEvals()
        if( outEvals ) {
            keys.add("eval_outputs")
            keys.add(computeEvalOutputCommands(outEvals))
        }

        // add variables referenced in the task script but not declared as input/output
        def vars = getTaskGlobalVars()
        if( vars ) {
            log.trace "Task: ${task.processor.name} > Adding script vars hash code: ${vars}"
            keys.add(vars.entrySet())
        }

        // add bin scripts referenced in the task script
        final binEntries = getTaskBinEntries(task.source)
        if( binEntries ) {
            log.trace "Task: ${task.processor.name} > Adding scripts on project bin path: ${-> binEntries.join('; ')}"
            keys.addAll(binEntries)
        }

        // add environment modules (`module` directive)
        final modules = task.getConfig().getModule()
        if( modules ) {
            keys.addAll(modules)
        }

        // add conda packages (`conda` directive)
        final conda = task.getCondaEnv()
        if( conda ) {
            keys.add(conda)
        }

        // add spack packages (`spack` and `arch` directives)
        final spack = task.getSpackEnv()
        final arch = task.getConfig().getArchitecture()

        if( spack ) {
            keys.add(spack)

            if( arch ) {
                keys.add(arch)
            }
        }

        // add stub run marker if enabled
        if( session.stubRun && task.config.getStubBlock() ) {
            keys.add('stub-run')
        }

        // compute task hash
        final mode = task.processor.getConfig().getHashMode()
        final hash = computeHash(keys, mode)

        // log task hash entries if enabled
        if( session.dumpHashes ) {
            session.dumpHashes == 'json'
                ? dumpHashEntriesJson(task, keys, mode, hash)
                : dumpHashEntriesLegacy(task, keys, mode, hash)
        }

        return hash
    }

    /**
     * Compute a deterministic string representation of eval output commands for cache hashing.
     */
    protected static String computeEvalOutputCommands(Map<String, String> outEvals) {
        assert outEvals != null && !outEvals.isEmpty(), "Eval outputs should not be null or empty"

        final result = new StringBuilder()
        final sortedEntries = outEvals.entrySet().sort { a, b -> a.key.compareTo(b.key) }

        for( final entry : sortedEntries ) {
            if( result.length() > 0 ) {
                result.append('\n')
            }
            result.append(entry.key).append('=').append(entry.value)
        }

        return result.toString()
    }

    protected Map<String,Object> getTaskGlobalVars() {
        return task.getTaskGlobalVars()
    }

    @Memoized
    protected List<Path> getTaskBinEntries(String script) {
        return task.getTaskBinEntries(script)
    }

    private String safeTaskName(TaskRun task) {
        return task != null ? task.lazyName() : task.processor.name
    }

    protected HashCode computeHash(List keys, CacheHelper.HashMode mode) {
        try {
            return createHashBuilder().withMode(mode).with(keys).build()
        }
        catch (Throwable e) {
            final msg = "Something went wrong while creating task hash for process '${task.processor.name}' -- Offending keys: ${ keys.collect { k -> "\n - type=${k.getClass().getName()} value=$k" } }"
            throw new UnexpectedException(msg,e)
        }
    }

    private void dumpHashEntriesJson(TaskRun task, List entries, CacheHelper.HashMode mode, hash) {
        final collector = (item) -> [
            hash: createHashBuilder().withMode(mode).with(item).build().toString(),
            type: item?.getClass()?.getName(),
            value: item?.toString()
        ]
        final json = JsonOutput.toJson(entries.collect(collector))
        log.info "[${safeTaskName(task)}] cache hash: ${hash}; mode: ${mode}; entries: ${JsonOutput.prettyPrint(json)}"
    }

    private void dumpHashEntriesLegacy(TaskRun task, List entries, CacheHelper.HashMode mode, hash) {
        final buffer = new StringBuilder()
        buffer.append("[${safeTaskName(task)}] cache hash: ${hash}; mode: $mode; entries: \n")
        for( final entry : entries ) {
            buffer.append( "  ${createHashBuilder().withMode(mode).with(entry).build()} [${entry?.getClass()?.getName()}] $entry \n")
        }
        log.info(buffer.toString())
    }
}
