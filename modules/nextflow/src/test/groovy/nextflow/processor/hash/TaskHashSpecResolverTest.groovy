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

import nextflow.SysEnv
import spock.lang.Specification

class TaskHashSpecResolverTest extends Specification {

    def cleanup() {
        SysEnv.pop()
    }

    def 'no version requested means no spec, so the run keeps the inherited hasher'() {
        given:
        SysEnv.push([:])

        expect:
        TaskHashSpecResolver.requestedSpec() == null
    }

    def 'a requested version resolves to its spec'() {
        given:
        SysEnv.push([NXF_TASK_HASH_VER: 'std/v2'])

        expect:
        TaskHashSpecResolver.requestedSpec().is(StdSpecs.STD_V2)
    }

    def 'an unknown version fails loudly rather than falling back'() {
        given:
        SysEnv.push([NXF_TASK_HASH_VER: 'std/nope'])

        when:
        TaskHashSpecResolver.requestedSpec()
        then:
        thrown(IllegalArgumentException)
    }

    def 'the factory abstains when no version is requested'() {
        given:
        SysEnv.push([:])

        expect:
        new SpecTaskHasherFactory().create(Mock(nextflow.processor.TaskRun)) == null
    }
}
