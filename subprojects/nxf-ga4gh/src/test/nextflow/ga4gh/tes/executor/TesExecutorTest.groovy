/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
