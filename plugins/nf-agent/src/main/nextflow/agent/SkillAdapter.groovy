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

import dev.langchain4j.skills.Skill
import dev.langchain4j.skills.Skills
import dev.langchain4j.skills.SkillResource as LcSkillResource
import groovy.transform.CompileStatic

/**
 * Maps core's portable, langchain4j-free {@link SkillDescriptor}s onto langchain4j
 * {@code Skill}s and bundles them into a {@code Skills} container (Tool Mode). The
 * resulting {@code Skills} provides both the tool provider (the {@code activate_skill}
 * / {@code read_skill_resource} tools) and the available-skills catalog the runner
 * injects into the system message.
 *
 * Mirrors {@link ModuleToolAdapter}: core does all filesystem work and passes plain
 * DTOs; the plugin builds the langchain4j objects, touching no filesystem path.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class SkillAdapter {

    /**
     * Build a langchain4j {@code Skills} container from the given portable descriptors.
     */
    static Skills toSkills(List<SkillDescriptor> descriptors) {
        final List<Skill> skills = new ArrayList<>()
        for( final SkillDescriptor d : descriptors ) {
            final List<LcSkillResource> resources = new ArrayList<>()
            if( d.resources ) {
                for( final SkillResource r : d.resources )
                    resources.add(LcSkillResource.builder().relativePath(r.relativePath).content(r.content).build())
            }
            skills.add(Skill.builder()
                .name(d.name)
                .description(d.description)
                .content(d.content)
                .resources(resources)
                .build())
        }
        return Skills.from(skills)
    }
}
