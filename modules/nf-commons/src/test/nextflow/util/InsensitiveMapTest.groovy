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

package nextflow.util

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class InsensitiveMapTest extends Specification {

    def 'should get value by case insensitive keys' () {
        given:
        def map = InsensitiveMap.of([alpha: 1, BETA: 2])

        expect:
        map.alpha == 1
        map.ALPHA == 1
        map.Alpha == 1
        and:
        map.get('alpha') == 1
        map.get('ALPHA') == 1
        map.get('Alpha') == 1
        and:
        map.beta == 2
        map.BETA == 2
        map.Beta == 2
        and:
        map.get('beta') == 2
        map.get('BETA') == 2
        map.get('Beta') == 2
        and:
        map.foo == null
        and:
        map.containsKey('alpha')
        map.containsKey('ALPHA')
        and:
        map.containsKey('beta')
        map.containsKey('BETA')
        and:
        !map.containsKey('foo')
    }

}
