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

package nextflow.extension.op

import groovy.transform.Canonical
import nextflow.prov.OperatorRun

/**
 * Associated a data value acquired by an operator with the corresponding
 * {@link OperatorRun} instance.
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
class OpDatum {
    Object value
    OperatorRun run

    static OpDatum of(Object value, OperatorRun run) {
        new OpDatum(value,run)
    }

    static Object unwrap(Object obj, List inputs=null) {
        if( obj instanceof Collection ) {
            return obj.collect(it-> unwrap(it,inputs))
        }
        if( obj instanceof OpDatum ) {
            if(inputs!=null)
                inputs.addAll(obj.run.inputIds)
            return obj.value
        }
        else
            return obj
    }
}
