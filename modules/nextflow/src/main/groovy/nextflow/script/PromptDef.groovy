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

    PromptDef(Closure closure, String source) {
        this.closure = closure
        this.source = source
    }

    @Override
    PromptDef clone() {
        (PromptDef) super.clone()
    }
}
