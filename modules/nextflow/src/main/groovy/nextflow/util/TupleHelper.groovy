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

package nextflow.util

import groovy.transform.CompileStatic

/**
 * SimpleHelper for tuple object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TupleHelper {

    static List listOf(x) {
        def result = new ArrayList(1)
        result.add(x)
        return result
    }

    static List listOf(x, y) {
        def result = new ArrayList(2)
        result.add(x)
        result.add(y)
        return result
    }

    static List listOf(x, y, z) {
        def result = new ArrayList(3)
        result.add(x)
        result.add(y)
        result.add(z)
        return result
    }

    static List listOf(x1, x2, x3, x4) {
        def result = new ArrayList(4)
        result.add(x1)
        result.add(x2)
        result.add(x3)
        result.add(x4)
        return result
    }

    static List listOf(x1, x2, x3, x4, x5) {
        def result = new ArrayList(5)
        result.add(x1)
        result.add(x2)
        result.add(x3)
        result.add(x4)
        result.add(x5)
        return result
    }

    static List listOf(x1, x2, x3, x4, x5, x6) {
        def result = new ArrayList(6)
        result.add(x1)
        result.add(x2)
        result.add(x3)
        result.add(x4)
        result.add(x5)
        result.add(x6)
        return result
    }

    static List listOf(x1, x2, x3, x4, x5, x6, x7) {
        def result = new ArrayList(7)
        result.add(x1)
        result.add(x2)
        result.add(x3)
        result.add(x4)
        result.add(x5)
        result.add(x6)
        result.add(x7)
        return result
    }

    static List listOf(x1, x2, x3, x4, x5, x6, x7, x8) {
        def result = new ArrayList(8)
        result.add(x1)
        result.add(x2)
        result.add(x3)
        result.add(x4)
        result.add(x5)
        result.add(x6)
        result.add(x7)
        result.add(x8)
        return result
    }
}
