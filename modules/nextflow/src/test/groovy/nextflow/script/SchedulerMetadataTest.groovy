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

import nextflow.Session
import nextflow.SysEnv
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class SchedulerMetadataTest extends Specification {

    @Unroll
    def 'should set enabled from the resolved executor' () {
        given:
        SysEnv.push(ENV)
        def session = Mock(Session) { getConfig()>>OPTS }

        expect:
        new SchedulerMetadata(session).enabled == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        OPTS                            | ENV                       | EXPECTED
        [:]                             | [:]                       | false
        [process:[executor:'local']]    | [:]                       | false
        [process:[executor:'seqera']]   | [:]                       | true
        [executor:[name:'seqera']]      | [:]                       | true
        [executor:[name:'awsbatch']]    | [:]                       | false
        [:]                             | [NXF_EXECUTOR:'seqera']   | true
        [:]                             | [NXF_EXECUTOR:'local']    | false
        // process.executor takes precedence over executor.name
        [process:[executor:'seqera'], executor:[name:'local']] | [:] | true
    }

    def 'should set enabled from the executor name' () {
        expect:
        new SchedulerMetadata('seqera').enabled
        and:
        !new SchedulerMetadata('local').enabled
        !new SchedulerMetadata((String)null).enabled
    }

    def 'should hold a mutable run id' () {
        given:
        def meta = new SchedulerMetadata('seqera')

        expect:
        meta.runId == null

        when:
        meta.runId = 'run-123'
        then:
        meta.runId == 'run-123'
    }

    def 'toMap should omit the run id until it is assigned' () {
        given:
        def meta = new SchedulerMetadata('seqera')

        expect: 'unset run id (e.g. at begin) -> only enabled'
        meta.toMap() == [enabled: true]

        when: 'run id assigned (e.g. by complete time)'
        meta.runId = 'run-123'
        then:
        meta.toMap() == [enabled: true, runId: 'run-123']
    }

}
