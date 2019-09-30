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
}
