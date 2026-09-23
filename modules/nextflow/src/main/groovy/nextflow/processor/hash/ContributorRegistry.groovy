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

import groovy.transform.CompileStatic

/**
 * Resolves a contributor by the canonical name a spec file refers to it by.
 *
 * A spec file describes composition — which keys, in what order, under which encoding —
 * but the extraction itself is code. This registry is the seam between the two, and the
 * reason a version that only adds, removes or reorders keys needs no Nextflow release
 * while one that needs a new extractor does.
 *
 * An unknown name is always an error. Falling back to a default would hash a task under
 * a spec nobody wrote.
 */
@CompileStatic
class ContributorRegistry {

    private static final Map<String,Contributor> REGISTRY = new LinkedHashMap<String,Contributor>()

    static {
        for( Contributor c : [
                Contributors.SESSION_ID,
                Contributors.PROCESS_NAME,
                Contributors.TASK_SOURCE,
                Contributors.CONTAINER,
                Contributors.INPUTS_RAW,
                Contributors.EVAL_OUTPUTS_RAW_MAP,
                Contributors.EVAL_OUTPUTS_DERIVED_STRING,
                Contributors.SCRIPT_VARS,
                Contributors.BIN_ENTRIES,
                Contributors.MODULE_BUNDLE,
                Contributors.ENV_MODULES,
                Contributors.CONDA,
                Contributors.SPACK_AND_ARCH,
                Contributors.STUB_MARKER ] ) {
            register(c)
        }
    }

    /**
     * Register a contributor under its canonical name. Plugins call this at startup to
     * make their own extractors referable from a spec file.
     */
    static void register(Contributor contributor) {
        final name = contributor.canonicalName()
        final existing = REGISTRY.get(name)
        if( existing != null && existing !== contributor ) {
            throw new IllegalArgumentException("Duplicate task hash contributor name: ${name}")
        }
        REGISTRY.put(name, contributor)
    }

    static Contributor get(String name) {
        final result = REGISTRY.get(name)
        if( result == null ) {
            throw new IllegalArgumentException("Unknown task hash contributor: ${name} -- available: ${REGISTRY.keySet().join(', ')}")
        }
        return result
    }

    static Set<String> names() {
        return Collections.unmodifiableSet(REGISTRY.keySet())
    }
}
