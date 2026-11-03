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

package nextflow.dataflow.ops

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.script.types.Bag
import nextflow.util.ArrayTuple
import nextflow.util.HashBag
/**
 * Implements the `groupBy` operator for typed workflows.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class GroupByOp {

    private DataflowReadChannel source

    GroupByOp(DataflowReadChannel source) {
        this.source = source
    }

    DataflowWriteChannel apply() {
        final target = CH.create()
        final groups = new HashMap<?,Bag<?>>()
        final sizes = new HashMap<?,Integer>()
        final emitted = new HashSet<>()

        final onNext = { ksv ->
            // obtain key, group size, and value
            def key, size, value
            if( ksv instanceof List && ksv.size() == 2 ) {
                (key, size, value) = [ ksv[0], -1, ksv[1] ]
            }
            else if( ksv instanceof List && ksv.size() == 3 ) {
                (key, size, value) = [ ksv[0], ksv[1] as Integer, ksv[2] ]
            }
            else {
                throw new ScriptRuntimeException("Operator `groupBy` expected a 3-tuple of (key, size, value) or a 2-tuple of (key, value) but received: ${ksv}")
            }

            if( emitted.contains(key) )
                throw new ScriptRuntimeException("Operator `groupBy` received too many values for grouping key: ${key} (expected ${size})")

            // append value to group
            final group = groups.computeIfAbsent(key, (k) -> new HashBag<>())
            group.add(value)

            // set group size
            if( sizes.containsKey(key) ) {
                if( size != sizes[key] )
                    throw new ScriptRuntimeException("Operator `groupBy` received inconsistent group size for key ${key} -- ${size} != ${sizes[key]}\n")
            }
            else {
                sizes[key] = size
            }

            // emit group when it is complete
            if( group.size() == size ) {
                target << new ArrayTuple<>([key, group])
                groups.remove(key)
                emitted.add(key)
            }
        }

        final onComplete = {
            // emit remaining groups
            groups.each { key, values ->
                if( sizes[key] != -1 )
                    throw new ScriptRuntimeException("Operator `groupBy` received too few values for grouping keys: ${key}")
                target << new ArrayTuple<>([key, values])
            }

            target << CH.stop()
        }

        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        return target
    }

}
