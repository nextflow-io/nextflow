/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.extension

import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Channel
import spock.lang.Timeout
import test.Dsl2Spec

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class OperatorDsl2Test extends Dsl2Spec {

    def 'should test unique' () {
        when:
        def result = dsl_eval("""
            Channel.of(1,2,3).unique()
        """) as DataflowReadChannel
        then:
        result.unwrap() == 1
        result.unwrap() == 2
        result.unwrap() == 3
        result.unwrap() == Channel.STOP
    }

    def 'should test unique with value' () {
        when:
        def result = dsl_eval("""
            Channel.value(1).unique()
        """) as DataflowReadChannel
        then:
        result.unwrap() == 1
    }

    def 'should test unique with collect' () {
        when:
        def result = dsl_eval("""
            Channel.of( 'a', 'b', 'c')
              .collect()
              .unique()
             .view()
            """) as DataflowReadChannel
       then:
       result.unwrap() == ['a','b','c']
    }

}
