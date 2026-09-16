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
import nextflow.util.EncodingRules

/**
 * The standard task hash specs, one per historical behaviour of master's hash.
 *
 * These are frozen. Changing a published spec changes the hash of every task already
 * recorded under it; a new behaviour gets a new id.
 */
@CompileStatic
class StdSpecs {

    /** Before #6679 "Record types" (2026-03-09): legacy encoding, no module bundle key. */
    static final TaskHashSpec STD_V1 = new TaskHashSpec('std/v1', [
        new KeyBinding(HashKey.SESSION_ID, Contributors.SESSION_ID),
        new KeyBinding(HashKey.PROCESS_NAME, Contributors.PROCESS_NAME),
        new KeyBinding(HashKey.TASK_SOURCE, Contributors.TASK_SOURCE),
        new KeyBinding(HashKey.CONTAINER, Contributors.CONTAINER),
        new KeyBinding(HashKey.INPUTS, Contributors.INPUTS_RAW),
        new KeyBinding(HashKey.EVAL_OUTPUTS, Contributors.EVAL_OUTPUTS_DERIVED_STRING),
        new KeyBinding(HashKey.SCRIPT_VARS, Contributors.SCRIPT_VARS),
        new KeyBinding(HashKey.BIN_ENTRIES, Contributors.BIN_ENTRIES),
        new KeyBinding(HashKey.ENV_MODULES, Contributors.ENV_MODULES),
        new KeyBinding(HashKey.CONDA, Contributors.CONDA),
        new KeyBinding(HashKey.SPACK, Contributors.SPACK_AND_ARCH),
        new KeyBinding(HashKey.STUB_MARKER, Contributors.STUB_MARKER)
    ], EncodingRules.LEGACY)

    /** #6679 (2026-03-09) to #6914 (2026-07-17): record-types encoding, still no module bundle. */
    static final TaskHashSpec STD_V2 = new TaskHashSpec('std/v2', STD_V1.bindings, EncodingRules.RECORD_TYPES)

    /** #6914 (2026-07-17) to #7575 (2026-09-03): module bundle key added. */
    static final TaskHashSpec STD_V3 = new TaskHashSpec('std/v3', [
        new KeyBinding(HashKey.SESSION_ID, Contributors.SESSION_ID),
        new KeyBinding(HashKey.PROCESS_NAME, Contributors.PROCESS_NAME),
        new KeyBinding(HashKey.TASK_SOURCE, Contributors.TASK_SOURCE),
        new KeyBinding(HashKey.CONTAINER, Contributors.CONTAINER),
        new KeyBinding(HashKey.INPUTS, Contributors.INPUTS_RAW),
        new KeyBinding(HashKey.EVAL_OUTPUTS, Contributors.EVAL_OUTPUTS_DERIVED_STRING),
        new KeyBinding(HashKey.SCRIPT_VARS, Contributors.SCRIPT_VARS),
        new KeyBinding(HashKey.BIN_ENTRIES, Contributors.BIN_ENTRIES),
        new KeyBinding(HashKey.MODULE_BUNDLE, Contributors.MODULE_BUNDLE),
        new KeyBinding(HashKey.ENV_MODULES, Contributors.ENV_MODULES),
        new KeyBinding(HashKey.CONDA, Contributors.CONDA),
        new KeyBinding(HashKey.SPACK, Contributors.SPACK_AND_ARCH),
        new KeyBinding(HashKey.STUB_MARKER, Contributors.STUB_MARKER)
    ], EncodingRules.RECORD_TYPES)

    /** Current master: post-#7575 (2026-09-03), eval hashed as a raw map. */
    static final TaskHashSpec STD_V4 = new TaskHashSpec('std/v4', [
        new KeyBinding(HashKey.SESSION_ID, Contributors.SESSION_ID),
        new KeyBinding(HashKey.PROCESS_NAME, Contributors.PROCESS_NAME),
        new KeyBinding(HashKey.TASK_SOURCE, Contributors.TASK_SOURCE),
        new KeyBinding(HashKey.CONTAINER, Contributors.CONTAINER),
        new KeyBinding(HashKey.INPUTS, Contributors.INPUTS_RAW),
        new KeyBinding(HashKey.EVAL_OUTPUTS, Contributors.EVAL_OUTPUTS_RAW_MAP),
        new KeyBinding(HashKey.SCRIPT_VARS, Contributors.SCRIPT_VARS),
        new KeyBinding(HashKey.BIN_ENTRIES, Contributors.BIN_ENTRIES),
        new KeyBinding(HashKey.MODULE_BUNDLE, Contributors.MODULE_BUNDLE),
        new KeyBinding(HashKey.ENV_MODULES, Contributors.ENV_MODULES),
        new KeyBinding(HashKey.CONDA, Contributors.CONDA),
        new KeyBinding(HashKey.SPACK, Contributors.SPACK_AND_ARCH),
        new KeyBinding(HashKey.STUB_MARKER, Contributors.STUB_MARKER)
    ], EncodingRules.RECORD_TYPES)

    static final TaskHashSpec DEFAULT = STD_V4

    static List<TaskHashSpec> all() {
        return [STD_V1, STD_V2, STD_V3, STD_V4]
    }

    static TaskHashSpec byId(String id) {
        final found = all().find { TaskHashSpec it -> it.id == id }
        if( !found ) {
            throw new IllegalArgumentException("Unknown task hash spec: ${id} -- available: ${all()*.id.join(', ')}")
        }
        return found
    }
}
