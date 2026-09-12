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

import groovy.transform.Canonical
import groovy.transform.CompileStatic

/**
 * Portable, langchain4j-free descriptor of a single file bundled with an agent
 * {@link SkillDescriptor skill} (e.g. a file under the skill's {@code references/}
 * directory). The content is loaded eagerly by core so the {@code nf-agent} plugin
 * can map it onto a langchain4j {@code SkillResource} without touching the
 * filesystem (the LLM reads it on demand via the {@code read_skill_resource} tool).
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
@CompileStatic
class SkillResource {
    /** Path of the resource relative to the skill directory (e.g. {@code references/guide.md}). */
    String relativePath
    /** The resource content, loaded eagerly by core. */
    String content
}
