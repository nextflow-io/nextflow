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
