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
package nextflow.script

import groovy.transform.CompileStatic

/**
 * Models the `prompt:` block of an agent definition. Mirrors {@link BodyDef}
 * but minimal: the prompt template is captured as a closure (evaluated per
 * invocation with the agent inputs in scope) plus its source text.
 */
@CompileStatic
class PromptDef implements Cloneable {

    final Closure closure
    final String source

    /**
     * The free-variable references (e.g. {@code params.*}, {@code task.ext.*})
     * captured from the prompt closure at parse time. These are folded into the
     * synthetic task {@link BodyDef#valRefs} so a prompt that closes over a
     * workflow-global enters the resume cache key (design §7.2/D3).
     */
    final List<TokenValRef> valRefs

    PromptDef(Closure closure, String source) {
        this(closure, source, Collections.<TokenValRef>emptyList())
    }

    PromptDef(Closure closure, String source, List<TokenValRef> valRefs) {
        this.closure = closure
        this.source = source
        this.valRefs = valRefs != null ? valRefs : Collections.<TokenValRef>emptyList()
    }

    /**
     * The names of the prompt's free-variable references (mirrors
     * {@link BodyDef#getValNames}).
     */
    List<String> getValNames() {
        valRefs*.name
    }

    @Override
    PromptDef clone() {
        (PromptDef) super.clone()
    }
}
