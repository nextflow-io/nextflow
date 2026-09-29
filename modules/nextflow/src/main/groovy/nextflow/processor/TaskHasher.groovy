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
package nextflow.processor

import java.nio.file.Path

import com.google.common.hash.HashCode
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.UnexpectedException
import nextflow.util.CacheHelper
/**
 * Implement task hash computation
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskHasher {

    private TaskRun task

    private TaskProcessor processor

    private Session session

    public TaskHasher(TaskRun task) {
        this.task = task
        this.processor = task.processor
        this.session = task.processor.session
    }

    public HashCode compute() {

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
        //
        // note: the map is hashed as an unordered collection of entries, so the
        // hash does not depend on the map iteration order, while still binding
        // each eval output name to its command
        final outEvals = task.getOutputEvals()
        if( outEvals ) {
            keys.add("eval_outputs")
            keys.add(outEvals)
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

        // add the fingerprint of the module resources bundle (e.g. module binaries)
        // since these files are not staged individually like the project bin dir
        final moduleBundle = session.enableModuleBinaries() ? processor.getModuleBundle() : null
        if( moduleBundle && moduleBundle.hasEntries() ) {
            final fingerprint = moduleBundle.fingerprint()
            log.trace "Task: ${task.processor.name} > Adding module bundle fingerprint: ${fingerprint}"
            keys.add(fingerprint)
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
     * Get the mapping of global variables that were referenced by
     * the task script, excluding references to `task.ext`.
     */
    Map<String,Object> getTaskGlobalVars() {
        final result = task.getGlobalVars(task.processor.getOwnerScript().getBinding())
        final directives = getTaskExtensionDirectiveVars()
        result.putAll(directives)
        return result
    }

    protected Map<String,Object> getTaskExtensionDirectiveVars() {
        final variableNames = task.getVariableNames()
        final result = new HashMap(variableNames.size())
        final taskConfig = task.config
        for( final key : variableNames ) {
            if( !key.startsWith('task.ext.') )
                continue
            final value = taskConfig.eval(key.substring(5))
            result.put(key, value)
        }

        return result
    }

    /**
     * This method scans the task command string looking for invocations of scripts
     * defined in the project bin folder.
     *
     * @param script The task command string
     * @return The list of paths of scripts in the project bin folder referenced in the task command
     */
    @Memoized
    List<Path> getTaskBinEntries(String script) {
        List<Path> result = []
        final tokenizer = new StringTokenizer(script, " \t\n\r\f()[]{};&|<>`")
        while( tokenizer.hasMoreTokens() ) {
            final token = tokenizer.nextToken()
            final path = session.binEntries.get(token)
            if( path )
                result.add(path)
        }
        return result
    }

    private String safeTaskName(TaskRun task) {
        return task != null ? task.lazyName() : task.processor.name
    }

    private HashCode computeHash(List keys, CacheHelper.HashMode mode) {
        try {
            return CacheHelper.hasher(keys, mode).hash()
        }
        catch (Throwable e) {
            final msg = "Something went wrong while creating task hash for process '${task.processor.name}' -- Offending keys: ${ keys.collect { k -> "\n - type=${k.getClass().getName()} value=$k" } }"
            throw new UnexpectedException(msg,e)
        }
    }

    private void dumpHashEntriesJson(TaskRun task, List entries, CacheHelper.HashMode mode, hash) {
        final collector = (item) -> [
            hash: CacheHelper.hasher(item, mode).hash().toString(),
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
            buffer.append( "  ${CacheHelper.hasher(entry, mode).hash()} [${entry?.getClass()?.getName()}] $entry \n")
        }
        log.info(buffer.toString())
    }
}
