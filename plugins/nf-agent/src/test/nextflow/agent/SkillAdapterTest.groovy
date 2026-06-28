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

import spock.lang.Specification

class SkillAdapterTest extends Specification {

    def 'should build a langchain4j Skills container with catalog and tool provider'() {
        given:
        def d = new SkillDescriptor('greet', 'a greeting skill', 'say hi politely',
            [new SkillResource('references/a.txt', 'AAA')])

        when:
        def skills = SkillAdapter.toSkills([d])

        then: 'a usable tool provider is produced'
        skills != null
        skills.toolProvider() != null

        and: 'the available-skills catalog advertises the name + description (what the model sees)'
        def catalog = skills.formatAvailableSkills()
        catalog.contains('greet')
        catalog.contains('a greeting skill')
    }

    def 'should handle a skill with no resources'() {
        when:
        def skills = SkillAdapter.toSkills([new SkillDescriptor('plain', 'no resources', 'body', null)])

        then:
        skills.formatAvailableSkills().contains('plain')
    }
}
