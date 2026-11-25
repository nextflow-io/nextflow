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
package nextflow.trace

import nextflow.Session
import nextflow.SysEnv
import spock.lang.Specification

class DefaultObserverFactoryTest extends Specification {

    private boolean hasResumeObserver(Map config) {
        def session = Stub(Session) {
            getConfig() >> config
        }
        return new DefaultObserverFactory().create(session).any { it instanceof ProgramResumeObserver }
    }

    def 'should not create the resume observer by default'() {
        given:
        SysEnv.push([:])

        expect:
        !hasResumeObserver([:])

        cleanup:
        SysEnv.pop()
    }

    def 'should create the resume observer when opted in by config'() {
        given:
        SysEnv.push([:])

        expect:
        hasResumeObserver([osc7577: true])

        cleanup:
        SysEnv.pop()
    }

    def 'should create the resume observer when opted in by the environment'() {
        given:
        SysEnv.push([NXF_OSC7577: 'true'])

        expect:
        hasResumeObserver([:])

        cleanup:
        SysEnv.pop()
    }

}
