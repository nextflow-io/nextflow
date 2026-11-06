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

import java.lang.reflect.Modifier

import groovy.json.JsonSlurper
import nextflow.lineage.model.v1beta1.TaskRun
import spock.lang.Specification
import spock.lang.Unroll

/**
 * Guards the key to lineage field mapping carried by the spec files.
 *
 * A consumer diffing two lineage records reads this mapping to know which fields belong
 * to which hash key. A key with no mapping drops out of the diff silently, and a mapping
 * naming a field that does not exist makes the consumer report no difference for a key it
 * never actually read -- both produce the dishonest answer the mapping exists to prevent.
 */
class SpecLineageFieldsTest extends Specification {

    static final Set<String> RECORD_FIELDS = TaskRun.declaredFields
            .findAll { !it.synthetic && !Modifier.isStatic(it.modifiers) }
            *.name as Set

    private static List<Map> keysOf(String id) {
        final path = "${StdSpecs.RESOURCE_DIR}/${id.replace('/', '-')}.json"
        final text = StdSpecs.getResourceAsStream(path).getText('UTF-8')
        return (new JsonSlurper().parseText(text) as Map).get('keys') as List<Map>
    }

    @Unroll
    def 'spec #id maps every key onto lineage fields that exist'() {
        given:
        final keys = keysOf(id)

        expect: 'no key is left unmapped -- an empty list is how a key says it is not recorded'
        keys.findAll { !it.containsKey('lineage') }*.key == []

        and:
        keys.collectMany { it.get('lineage') as List<String> }.toSet() - RECORD_FIELDS == [] as Set

        where:
        id << StdSpecs.IDS
    }
}
