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

package nextflow.cloud.azure.batch

import nextflow.processor.TaskProcessor
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzJobKeyTest extends Specification {

    def 'should validate equals and hashcode' () {
        given:
        def p1 = Mock(TaskProcessor)
        def p2 = Mock(TaskProcessor)
        def k1 = new AzJobKey(p1, 'foo')
        def k2 = new AzJobKey(p1, 'foo')
        def k3 = new AzJobKey(p2, 'foo')
        def k4 = new AzJobKey(p1, 'bar')

        expect:
        k1 == k2
        k1 != k3
        k1 != k4
        and:
        k1.hashCode() == k2.hashCode()
        k1.hashCode() != k3.hashCode()
        k1.hashCode() != k4.hashCode()
    }

}
