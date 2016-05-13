/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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

package nextflow.extension

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.Session
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowMergeExtensionTest extends Specification {

    @Shared
    Session session

    def setup() {
        session = new Session()
    }

    def cleanup() {
        assert !session.dag.isEmpty()
    }

    def 'should merge with open array'() {

        when:
        def alpha = Channel.from(1, 3, 5);
        def beta = Channel.from(2, 4, 6);
        def delta = Channel.from(7,8,1);

        def result = alpha.merge( beta, delta ) { a,b,c -> [a,b,c] }

        then:
        result instanceof DataflowQueue
        result.val == [1,2,7]
        result.val == [3,4,8]
        result.val == [5,6,1]
        result.val == Channel.STOP
    }

    def 'should merge with list'() {

        when:
        def alpha = Channel.from(1, 3, 5);
        def beta = Channel.from(2, 4, 6);
        def delta = Channel.from(7,8,1);

        def result = alpha.merge( [beta, delta] ) { a,b,c -> [c,b,a] }

        then:
        result instanceof DataflowQueue
        result.val == [7,2,1]
        result.val == [8,4,3]
        result.val == [1,6,5]
        result.val == Channel.STOP
    }

    def 'should merge with queue'() {

        when:
        def alpha = Channel.from(1, 3, 5);
        def beta = Channel.from(2, 4, 6);

        def result = alpha.merge(beta) { a,b -> [a, b+1] }

        then:
        result instanceof DataflowQueue
        result.val == [1,3]
        result.val == [3,5]
        result.val == [5,7]
        result.val == Channel.STOP
    }

    def 'should merge with variables'() {

        when:
        def alpha = Channel.value('Hello');
        def beta = Channel.value('World')

        def result = alpha.merge(beta) { a,b -> [a, b] }

        then:
        result instanceof DataflowVariable
        result.val == ['Hello','World']
        result.val == ['Hello','World']
    }

}
