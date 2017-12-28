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

package nextflow.extension

import nextflow.Session
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowHelperTest extends Specification {


    def setupSpec() {
        new Session()
    }

    def 'should subscribe handlers'() {

        when:
        DataflowHelper.checkSubscribeHandlers( [:] )
        then:
        thrown(IllegalArgumentException)

        when:
        DataflowHelper.checkSubscribeHandlers( [ onNext:{}] )
        then:
        true

        when:
        DataflowHelper.checkSubscribeHandlers( [ onNext:{}, xxx:{}] )
        then:
        thrown(IllegalArgumentException)

        when:
        DataflowHelper.checkSubscribeHandlers( [ xxx:{}] )
        then:
        thrown(IllegalArgumentException)
    }

    @Unroll
    def 'should split entry' () {
        when:
        def pair = DataflowHelper.split(pivot, entry)
        then:
        pair.keys == keys
        pair.values == values

        where:
        pivot           | entry                         | keys          | values
        [0]             | ['A','B','C','D','E','F']     | ['A']         | ['B','C','D','E','F']
        [0,1]           | ['A','B','C','D','E','F']     | ['A','B']     | ['C','D','E','F']
        [0,2]           | ['A','B','C','D','E','F']     | ['A','C']     | ['B','D','E','F']
        [0,1,4]         | ['A','B','C','D','E','F']     | ['A','B','E'] | ['C','D','F']
        [0]             | 'A'                           | ['A']         | []
    }
}
