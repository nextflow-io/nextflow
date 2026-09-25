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

import com.google.common.hash.HashCode
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.UnexpectedException
import nextflow.processor.TaskHasher
import nextflow.processor.TaskRun
import nextflow.util.CacheHelper
import nextflow.util.HashBuilder

/**
 * A {@link TaskHasher} that computes the task hash from a {@link TaskHashSpec} — an
 * ordered list of named key bindings plus encoding rules — rather than from a hard-coded
 * key sequence.
 *
 * It is opt-in: {@link SpecTaskHasherFactory} returns one only when a task hash version
 * is configured, so a run that asks for nothing keeps the inherited {@link TaskHasher}
 * behaviour byte for byte. There are no subclasses of this class — a different hash
 * version, or a plugin's own hashing, is a different spec value, not a new type.
 */
@Slf4j
@CompileStatic
class SpecTaskHasher extends TaskHasher {

    private final HashContext ctx

    private final TaskHashSpec spec

    SpecTaskHasher(HashContext ctx, TaskHashSpec spec) {
        super(ctx.task)
        this.ctx = ctx
        this.spec = spec
    }

    SpecTaskHasher(TaskRun task, TaskHashSpec spec) {
        this(new HashContext(task), spec)
    }

    TaskHashSpec getSpec() {
        return spec
    }

    /** The flat, ordered value list the spec produces for this task. */
    List<Object> collectKeys() {
        final keys = new ArrayList<Object>()
        for( KeyBinding binding : spec.bindings ) {
            keys.addAll(binding.contributor.emit(ctx))
        }
        return keys
    }

    @Override
    HashCode compute() {
        final keys = collectKeys()
        final mode = ctx.task.processor.getConfig().getHashMode()
        try {
            return spec.encoding
                .apply(new HashBuilder().withHasher(HashBuilder.defaultHasher()).withMode(mode))
                .with(keys)
                .build()
        }
        catch( Throwable e ) {
            final msg = "Something went wrong while creating task hash for process '${ctx.processor.name}' under spec '${spec.id}' -- Offending keys: ${ keys.collect { k -> "\n - type=${k?.getClass()?.getName()} value=$k" } }"
            throw new UnexpectedException(msg, e)
        }
    }

    /**
     * Per-key digest of this task under this spec, for explaining a cache miss.
     *
     * Strictly additive: it recomputes nothing that feeds compute(), and a key that
     * emits no value is omitted rather than digested as empty.
     */
    @CompileStatic
    Map<HashKey,HashCode> explain() {
        final mode = ctx.task.processor.getConfig().getHashMode()
        final result = new LinkedHashMap<HashKey,HashCode>()
        for( KeyBinding binding : spec.bindings ) {
            final values = binding.contributor.emit(ctx)
            if( !values ) {
                continue
            }
            result.put(binding.key, spec.encoding
                .apply(new HashBuilder().withHasher(HashBuilder.defaultHasher()).withMode(mode))
                .with(values)
                .build())
        }
        return result
    }

    /**
     * Named per-key entries for `-dump-hashes json`, prefixed by the spec identity.
     */
    String dumpJson() {
        final entries = new ArrayList<Map<String,Object>>()
        entries.add([spec: spec.id, fingerprint: spec.fingerprint()] as Map<String,Object>)
        for( Map.Entry<HashKey,HashCode> e : explain().entrySet() ) {
            entries.add([key: e.key.name(), hash: e.value.toString()] as Map<String,Object>)
        }
        return JsonOutput.prettyPrint(JsonOutput.toJson(entries))
    }

    /**
     * Per-entry dump for `-dump-hashes` (legacy format), iterating the flat key list.
     */
    String dumpLegacy(HashCode hash) {
        final mode = ctx.task.processor.getConfig().getHashMode()
        final keys = collectKeys()
        final buffer = new StringBuilder()
        buffer.append("[${ctx.task.lazyName()}] cache hash: ${hash}; mode: $mode; entries: \n")
        for( final entry : keys ) {
            buffer.append("  ${CacheHelper.hasher(entry, mode).hash()} [${entry?.getClass()?.getName()}] $entry \n")
        }
        return buffer.toString()
    }
}
