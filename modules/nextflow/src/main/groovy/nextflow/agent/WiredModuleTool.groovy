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
package nextflow.agent

import groovy.transform.CompileStatic
import groovy.transform.TupleConstructor
import io.seqera.npr.api.schema.v1.ModuleMetadata
import nextflow.module.ModuleSpec
import nextflow.script.ProcessDef

/**
 * Everything known about ONE brokered {@code nf:module_run} tool, gathered while the agent is
 * lowered by {@link ModuleToolResolver} and handed to {@link ModuleToolBridge} to wire.
 *
 * <p>The facts travel together rather than in parallel name-keyed maps because they are only ever
 * read together: the bridge asks "does this tool have a spec, and is there registry metadata to
 * describe it with" once per tool, in list order.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@TupleConstructor
@CompileStatic
class WiredModuleTool {

    /** The wire name the model sees, i.e. the process name selected through {@code nf:module_run}. */
    final String name

    /** The process template cloned for every tool call. */
    final ProcessDef proc

    /**
     * The sibling {@code meta.yml} spec, or {@code null} for a locally-defined process. Its
     * presence is what selects spec-driven marshalling over the scalar typed-I/O path.
     */
    final ModuleSpec spec

    /**
     * The public registry {@link ModuleMetadata}, or {@code null} when the module is not a
     * registry install or its metadata could not be fetched. When present it is the descriptor
     * source (description + input schema); marshalling stays on the spec + ProcessDef.
     */
    final ModuleMetadata metadata

    /** Whether the metadata's module is nf-core scoped (for the {@code meta.id} convention);
     * always {@code false} when there is no {@link #metadata}. */
    final boolean nfCore

}
