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
 * Portable, langchain4j-free descriptor of an Anthropic-style agent <em>skill</em>
 * (a {@code SKILL.md} folder). Mirrors {@link ToolDescriptor}: core resolves and
 * parses the skill (locally or after a remote fetch) into this plain DTO, and the
 * {@code nf-agent} plugin maps it onto a langchain4j {@code Skill} — so core stays
 * free of any LLM client dependency and the plugin never touches a filesystem path.
 *
 * <p>{@code name}/{@code description} come from the {@code SKILL.md} YAML frontmatter
 * (the model sees them in the available-skills catalog); {@code content} is the
 * markdown body returned by the {@code activate_skill} tool; {@code resources} are
 * the bundled files returned by {@code read_skill_resource}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
@CompileStatic
class SkillDescriptor {
    String name
    String description
    String content
    List<SkillResource> resources
}
