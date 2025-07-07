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

package nextflow.container.resolver

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ContainerInfoTest extends Specification {

    def 'should create container info' () {
        when:
        def info = new ContainerInfo('a','b','c')
        then:
        info.source == 'a'
        info.target == 'b'
        info.hashKey == 'c'
    }

    def 'should check truth' () {
        expect:
        new ContainerInfo('a')
        new ContainerInfo(null,'b')
        new ContainerInfo(null,null,'c')
        new ContainerInfo('a','b','c')
        and:
        !ContainerInfo.EMPTY
        !new ContainerInfo()
    }
}
