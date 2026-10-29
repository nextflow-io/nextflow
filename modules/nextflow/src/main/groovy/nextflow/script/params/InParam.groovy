/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.script.params

import groovyx.gpars.dataflow.DataflowReadChannel

/**
 * Basic interface for *all* input parameters
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface InParam extends Cloneable {

    String getName()

    DataflowReadChannel getInChannel()

    Object getRawChannel()

    // NOTE: declared as explicit getters (rather than `short index` /
    // `short mapIndex`) because Groovy 5 compiles a bare typed field in an
    // interface as a `static final` constant (value 0) instead of an abstract
    // property. That constant shadowed the real instance value during dynamic
    // property access (e.g. in `BaseInParam.decodeInputs`), breaking input
    // resolution.
    short getIndex()

    short getMapIndex()

    def decodeInputs( List values )

}
