/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.ga4gh.tes.executor

import nextflow.Session
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TesExecutorTest extends Specification {

    def 'should get endpoint' () {
        given:
        def session = Mock(Session)
        def exec = new TesExecutor(session: session)

        when:
        def result = exec.getEndPoint()
        then:
        session.getConfigAttribute('executor.tes.endpoint','http://localhost:8000') >> 'http://foo.com'
        result == 'http://foo.com'
    }


    def 'should resolve endpoint from env'() {
        given:
        def ENV = [NXF_EXECUTOR_TES_ENDPOINT: 'http://back.end.com' ]
        def session = Spy(Session)
        def exec = Spy(TesExecutor)

        when:
        def result = exec.getEndPoint()
        then:
        1 * exec.getSession() >> session
        1 * session.getSystemEnv() >> ENV
        result == 'http://back.end.com'

    }

}
