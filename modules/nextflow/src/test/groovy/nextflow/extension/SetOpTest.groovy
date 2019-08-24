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

import spock.lang.Timeout

import test.Dsl2Spec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class SetOpTest extends Dsl2Spec {

    def 'should set a channel in the global context' () {
        when:
        def result = dsl_eval(/
            Channel.from(1,2,3) | set { foo }
            foo | map { it *2 }
        /)
        then:
        result.val == 2
        result.val == 4
        result.val == 6

        when:
        result = dsl_eval(/
            Channel.value(5) | set { foo }
                foo | map { it *2 }
        /)
        then:
        result.val == 10
    }

    def 'should invoke set with dot notation' () {
        when:
        def result = dsl_eval(/
            Channel.from(1,2,3).set { foo } 
            foo.map { it *2 }
        /)
        then:
        result.val == 2
        result.val == 4
        result.val == 6

        when:
        result = dsl_eval(/
            Channel.value('hello').set { foo } 
            foo.map { it.toUpperCase() }
        /)
        then:
        result.val == 'HELLO'
    }


    def 'should assign multiple channels in the current binding' () {
        when:
        def result = dsl_eval(/
        def ch1 = Channel.value('X')
        def ch2 = Channel.value('Y')

        new nextflow.script.ChannelOut([ch1]) .set { foo }   
        return foo
        /)
        then:
        result.val == 'X'

        when:
        result = dsl_eval(/
        def ch1 = Channel.value('X')
        def ch2 = Channel.value('Y')

        new nextflow.script.ChannelOut([ch1, ch2]) .set { bar }
        return bar
        /)
        then:
        result[0].val == 'X'
        result[1].val == 'Y'
    }

}
