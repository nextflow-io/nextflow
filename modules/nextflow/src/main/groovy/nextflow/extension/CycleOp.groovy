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

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel

/**
 * Implements the cycle operator that pairs channel items with cycling indices
 *
 * @author [Your Name]
 */
@CompileStatic
class CycleOp {

    private DataflowReadChannel source
    private List indices
    private int currentIndex = 0

    CycleOp(DataflowReadChannel source, range) {
        this.source = source
        this.indices = range instanceof Range ? range.collect() : range as List
    }

    DataflowWriteChannel apply() {
        final target = CH.create()
        
        final events = [
            onNext: { item ->
                target.bind([indices[currentIndex], item])
                currentIndex = (currentIndex + 1) % indices.size()
            },
            onComplete: { target << Channel.STOP }
        ]
        
        DataflowHelper.subscribeImpl(source, events)
        return target
    }
} 