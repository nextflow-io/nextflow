/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.file

import spock.lang.Specification

import nextflow.Session
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PathVisitorTest extends Specification {

    def 'should return the same executor instance' () {

        given:
        def session1 = Mock(Session)
        def session2 = Mock(Session)

        when:
        def X = PathVisitor.createExecutor(session1)
        then:
        X != null

        when:
        def Y = PathVisitor.createExecutor(session1)
        then:
        X.is(Y)

        when:
        def P = PathVisitor.createExecutor(session2)
        then:
        !P.is(X)

        when:
        def Q = PathVisitor.createExecutor(session1)
        then:
        Q.is(X)
    }

}
