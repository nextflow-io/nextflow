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
 */

package nextflow.extension

import nextflow.Channel
import nextflow.Session
import spock.lang.Specification

/**
 * Test GroupKey preservation in join operations
 *
 * @author Your Name
 */
class JoinOpGroupKeyTest extends Specification {

    def setup() {
        new Session()
    }

    def 'should preserve GroupKey when joining channels' () {
        given:
        def key1 = new GroupKey('X', 2)
        def key2 = new GroupKey('Y', 3)
        
        def ch1 = Channel.of([key1, 1], [key2, 2])
        def ch2 = Channel.of(['X', 'a'], ['Y', 'b'])

        when:
        def op = new JoinOp(ch1, ch2)
        def result = op.apply().toList().getVal()
        
        then:
        result.size() == 2
        
        // Check that GroupKey is preserved in the output
        result.each { tuple ->
            assert tuple[0] instanceof GroupKey
            assert tuple.size() == 3
        }
        
        // Verify the actual values
        def sorted = result.sort { it[0].toString() }
        sorted[0][0].toString() == 'X'
        sorted[0][0].groupSize == 2
        sorted[0][1] == 1
        sorted[0][2] == 'a'
        
        sorted[1][0].toString() == 'Y'
        sorted[1][0].groupSize == 3
        sorted[1][1] == 2
        sorted[1][2] == 'b'
    }

    def 'should preserve GroupKey when GroupKey is on right channel' () {
        given:
        def key1 = new GroupKey('X', 2)
        def key2 = new GroupKey('Y', 3)
        
        def ch1 = Channel.of(['X', 'a'], ['Y', 'b'])
        def ch2 = Channel.of([key1, 1], [key2, 2])

        when:
        def op = new JoinOp(ch1, ch2)
        def result = op.apply().toList().getVal()
        
        then:
        result.size() == 2
        
        // Check that GroupKey is preserved in the output
        result.each { tuple ->
            assert tuple[0] instanceof GroupKey
            assert tuple.size() == 3
        }
    }

    def 'should handle mix of GroupKey and plain keys correctly' () {
        given:
        def key1 = new GroupKey('X', 2)
        
        def ch1 = Channel.of([key1, 1], ['Y', 2])  // Mix of GroupKey and plain key
        def ch2 = Channel.of(['X', 'a'], ['Y', 'b'])

        when:
        def op = new JoinOp(ch1, ch2)
        def result = op.apply().toList().getVal().sort { it[0].toString() }
        
        then:
        result.size() == 2
        
        // First tuple should have GroupKey
        result[0][0] instanceof GroupKey
        result[0][0].toString() == 'X'
        result[0][0].groupSize == 2
        
        // Second tuple should have plain string
        result[1][0] == 'Y'
        !(result[1][0] instanceof GroupKey)
    }
} 