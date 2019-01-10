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

package nextflow.k8s.model

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PodNodeSelectorTest extends Specification {

    def 'should create node selector' () {

        expect: 
        new PodNodeSelector(selector).toSpec() == spec

        where:
        selector            | spec
        ''                  | [:]
        'foo=1'             | [foo:'1']
        'x=a,y=2,z=9'       | [x:'a',y:'2',z:'9']
        'x= a , y=2 , z =9' | [x:'a',y:'2',z:'9']
        'gpu,intel'         | [gpu:'true',intel: 'true']
        [foo:1, bar: 'two'] | [foo:'1', bar:'two']
    }


}
