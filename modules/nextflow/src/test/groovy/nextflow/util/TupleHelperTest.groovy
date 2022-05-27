/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.util

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TupleHelperTest extends Specification {

    def 'check single value' () {
        expect:
        TupleHelper.listOf('a') == ['a']
        TupleHelper.listOf('a','b') == ['a', 'b']
        TupleHelper.listOf('a','b','z') == ['a', 'b', 'z']
        TupleHelper.listOf('1','2','3', '4') == '1'..'4'
        TupleHelper.listOf('1','2','3', '4', '5') == '1'..'5'
        TupleHelper.listOf('1','2','3', '4', '5', '6') == '1'..'6'
        TupleHelper.listOf('1','2','3', '4', '5', '6', '7') == '1'..'7'
        TupleHelper.listOf('1','2','3', '4', '5', '6', '7', '8') == '1'..'8'
    }

}
