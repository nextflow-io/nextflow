/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.extension

import static nextflow.extension.DataflowHelper.*

import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel

/**
 * Implements Nextflow Mix operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class MixOp {

    private DataflowReadChannel source
    private DataflowReadChannel[] others

    MixOp(DataflowReadChannel source, DataflowReadChannel other) {
        this.source = source
        this.others = [other]
    }

    MixOp(DataflowReadChannel source, DataflowReadChannel[] others) {
        this.source = source
        this.others = others
    }

    DataflowWriteChannel apply() {
        def target = CH.create()
        def count = new AtomicInteger( others.size()+1 )
        def handlers = [
                onNext: { target << it },
                onComplete: { if(count.decrementAndGet()==0) { target << Channel.STOP } }
        ]

        subscribeImpl(source, handlers)
        for( def it : others ) {
            subscribeImpl(it, handlers)
        }

        def allSources = [source]
        allSources.addAll(others)
        return target
    }

}
