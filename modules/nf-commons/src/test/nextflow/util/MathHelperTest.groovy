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
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MathHelperTest extends Specification {

    @Unroll
    def 'should validate ceil div int' () {
        expect:
        MathHelper.ceilDiv(X, Y) == Z
        where:
        X   | Y     | Z
        0   | 10    | 0
        1   | 10    | 1
        9   | 10    | 1
        10  | 10    | 1
        11  | 10    | 2
        19  | 10    | 2
        20  | 10    | 2
        21  | 10    | 3
    }

    @Unroll
    def 'should validate ceil div long' () {
        expect:
        MathHelper.ceilDiv(X, Y) == Z
        where:
        X                   | Y                 | Z
        0                   | 10_000_000_000    | 0
        1                   | 10_000_000_000    | 1
        10_000_000_000-1    | 10_000_000_000    | 1
        10_000_000_000      | 10_000_000_000    | 1
        and:
        10_000_000_001      | 10_000_000_000    | 2
        20_000_000_000      | 10_000_000_000    | 2
        and:
        20_000_000_001      | 10_000_000_000    | 3
        30_000_000_000      | 10_000_000_000    | 3
    }
}
