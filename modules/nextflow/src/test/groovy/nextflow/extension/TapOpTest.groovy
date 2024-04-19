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

import test.Dsl2Spec

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TapOpTest extends Dsl2Spec {

    def 'should forward a queue channel into multiple targets' () {
        when:
        def result = dsl_eval('''
            Channel.of(1, 2, 3).tap { foo ; bar }
            [ foo.collect(), bar.collect() ]
            ''')

        then:
        result[0].getVal() == [1, 2, 3]
        result[1].getVal() == [1, 2, 3]
    }

    def 'should forward a value channel into multiple targets' () {
        when:
        def result = dsl_eval('''
            Channel.value(4).tap { foo ; bar }
            [ foo, bar ]
            ''')
        then:
        result[0].getVal() == 4
        result[1].getVal() == 4
    }

}
