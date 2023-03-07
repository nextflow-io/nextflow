/*
 * Copyright 2023, Seqera Labs
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

import nextflow.Channel
import nextflow.Session
import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class GroupMapOpTest extends Specification {

    def setup() {
        new Session()
    }

    def 'should group items using dynamic group size' () {
        given:
        def k1 = new GroupKey('k1', 2)
        def k2 = new GroupKey('k2', 3)
        def k3 = new GroupKey('k3', 4)

        def items = [
            [group: k1, name: 'a'],
            [group: k1, name: 'b'],
            [group: k2, name: 'x'],
            [group: k3, name: 'q'],
            [group: k1, name: 'd'],
            [group: k1, name: 'c'],
            [group: k2, name: 'y'],
            [group: k1, name: 'f'],
            [group: k2, name: 'z']
        ]

        when:
        // here the size is defined as operator argument
        def result = items.channel().groupMap(by: 'group', size: 2)
        then:
        result.val == [group: k1, name: ['a', 'b']]
        result.val == [group: k1, name: ['d', 'c']]
        result.val == [group: k2, name: ['x', 'y']]
        result.val == Channel.STOP

        when:
        // here the size is inferred by the key itself
        result = items.channel().groupMap(by: 'group')
        then:
        result.val == [group: k1, name: ['a', 'b']]
        result.val == [group: k1, name: ['d', 'c']]
        result.val == [group: k2, name: ['x', 'y', 'z']]
        result.val == Channel.STOP

        when:
        result = items.channel().groupMap(by: 'group', remainder: true)
        then:
        result.val == [group: k1, name: ['a', 'b']]
        result.val == [group: k1, name: ['d', 'c']]
        result.val == [group: k2, name: ['x', 'y', 'z']]
        result.val == [group: k3, name: ['q']]
        result.val == [group: k1, name: ['f']]
        result.val == Channel.STOP
    }
}
