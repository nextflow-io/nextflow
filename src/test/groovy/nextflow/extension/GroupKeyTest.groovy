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

package nextflow.extension

import spock.lang.Specification

import nextflow.extension.GroupKey
import nextflow.util.KryoHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GroupKeyTest extends Specification {

    def 'should mimic target object' () {
        given:
        def hello = 'Hello world!'
        when:
        def group = new GroupKey(hello, 10)

        then:
        group.size() == hello.size()
        group == hello
        group.equals(hello)
        group.hashCode() == hello.hashCode()
        group.groupSize == 10

    }

    def 'should verify equals and hashCode' () {
        given:
        def hello = 'Hello world!'
        when:
        def g1 = new GroupKey(hello, 10)
        def g2 = new GroupKey(hello, 10)
        def g3 = new GroupKey(hello.reverse(), 10)

        then:
        g1 == g2
        g1.equals(g2)
        
        g1 != g3
        !g1.equals(g3)

        g1.hashCode() == g2.hashCode()
        g1.hashCode() != g3.hashCode()
    }

    def 'should serialise/de-serialise a groupkey' () {

        given:
        def hello = 'Hello world!'
        def group = new GroupKey(hello, 10)
        when:
        def buff = KryoHelper.serialize(group)
        then:
        def copy = (GroupKey)KryoHelper.deserialize(buff)
        and:
        copy == group

    }

}
