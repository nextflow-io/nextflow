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

import groovy.transform.CompileStatic
import nextflow.script.bundle.ResourcesBundle

/**
 * The contributor implementations, one per key site of TaskHasher.compute().
 *
 * Every body here mirrors the corresponding lines of that method exactly, including
 * how many values it appends and the fact that it appends none when the key is
 * absent. Byte-exactness lives in this file.
 */
@CompileStatic
class Contributors {

    static Contributor of(String name, Closure<List<Object>> fn) {
        return new Contributor() {
            @Override
            String canonicalName() {
                return name
            }

            @Override
            List<Object> emit(HashContext ctx) {
                return (List<Object>) fn.call(ctx)
            }
        }
    }

    static final Contributor SESSION_ID = of('sessionId') { HashContext ctx ->
        return [ctx.session.uniqueId] as List<Object>
    }

    static final Contributor PROCESS_NAME = of('processName') { HashContext ctx ->
        return [ctx.processor.name] as List<Object>
    }

    static final Contributor TASK_SOURCE = of('taskSource') { HashContext ctx ->
        return [ctx.task.source] as List<Object>
    }

    static final Contributor CONTAINER = of('containerFingerprint') { HashContext ctx ->
        if( !ctx.task.isContainerEnabled() ) {
            return [] as List<Object>
        }
        return [ctx.task.getContainerFingerprint()] as List<Object>
    }

    static final Contributor INPUTS_RAW = of('inputs.raw') { HashContext ctx ->
        final out = new ArrayList<Object>()
        for( final entry : ctx.task.inputs ) {
            out.add(entry.key.name)
            out.add(entry.value)
        }
        return out
    }

    /** Post-#7575: the eval map is hashed directly. */
    static final Contributor EVAL_OUTPUTS_RAW_MAP = of('evalOutputs.rawMap') { HashContext ctx ->
        final outEvals = ctx.task.getOutputEvals()
        if( !outEvals ) {
            return [] as List<Object>
        }
        return ['eval_outputs', outEvals] as List<Object>
    }

    /** Pre-#7575: a sorted "name=command" string, one entry per line. */
    static final Contributor EVAL_OUTPUTS_DERIVED_STRING = of('evalOutputs.derivedString') { HashContext ctx ->
        final outEvals = ctx.task.getOutputEvals()
        if( !outEvals ) {
            return [] as List<Object>
        }
        return ['eval_outputs', derivedEvalCommands(outEvals)] as List<Object>
    }

    static final Contributor SCRIPT_VARS = of('scriptVars') { HashContext ctx ->
        final vars = ctx.globalVars()
        if( !vars ) {
            return [] as List<Object>
        }
        return [vars.entrySet()] as List<Object>
    }

    static final Contributor BIN_ENTRIES = of('binEntries') { HashContext ctx ->
        final entries = ctx.binEntries()
        if( !entries ) {
            return [] as List<Object>
        }
        return new ArrayList<Object>(entries)
    }

    static final Contributor MODULE_BUNDLE = of('moduleBundleFingerprint') { HashContext ctx ->
        final ResourcesBundle bundle = ctx.session.enableModuleBinaries()
            ? ctx.processor.getModuleBundle()
            : null
        if( !bundle || !bundle.hasEntries() ) {
            return [] as List<Object>
        }
        return [bundle.fingerprint()] as List<Object>
    }

    static final Contributor ENV_MODULES = of('envModules') { HashContext ctx ->
        final modules = ctx.task.getConfig().getModule()
        if( !modules ) {
            return [] as List<Object>
        }
        return new ArrayList<Object>(modules)
    }

    static final Contributor CONDA = of('condaEnv') { HashContext ctx ->
        final conda = ctx.task.getCondaEnv()
        if( !conda ) {
            return [] as List<Object>
        }
        return [conda] as List<Object>
    }

    /** arch contributes only when spack is set — it is nested inside that branch today. */
    static final Contributor SPACK_AND_ARCH = of('spackEnvAndArch') { HashContext ctx ->
        final spack = ctx.task.getSpackEnv()
        if( !spack ) {
            return [] as List<Object>
        }
        final arch = ctx.task.getConfig().getArchitecture()
        if( !arch ) {
            return [spack] as List<Object>
        }
        return [spack, arch] as List<Object>
    }

    /**
     * TaskConfig.getStubBlock() is protected and lives in a different package
     * (nextflow.processor), so it is not reachable from here under @CompileStatic.
     * TaskRun.getStubSource() is the existing public equivalent: it is null iff
     * getStubBlock() is null, since a real stub closure's source is never null.
     */
    static final Contributor STUB_MARKER = of('stubMarker') { HashContext ctx ->
        if( ctx.session.stubRun && ctx.task.getStubSource() != null ) {
            return ['stub-run'] as List<Object>
        }
        return [] as List<Object>
    }

    /**
     * Reproduces TaskHasher.computeEvalOutputCommands(), removed by #7575.
     * Retained verbatim because specs std/v1..v3 depend on its exact output.
     */
    protected static String derivedEvalCommands(Map<String,String> outEvals) {
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
}
