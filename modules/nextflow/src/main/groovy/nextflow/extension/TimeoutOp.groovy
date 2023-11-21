/*
 * Copyright 2013-2023, Seqera Labs
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

import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.CopyChannelsClosure
import nextflow.Channel
import nextflow.extension.DataflowHelper
import nextflow.util.Duration

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TimeoutOp {

    private DataflowReadChannel source
    private Duration duration

    TimeoutOp(final DataflowReadChannel source, Duration duration) {
        this.source = source
        this.duration = duration
    }

    DataflowWriteChannel apply() {
        final target = CH.create()
        final op = DataflowHelper.newOperator(
                source,
                target,
                new ChainWithClosure(new CopyChannelsClosure()))

        final task = {
            op.bindOutput(Channel.STOP)
            op.terminate()
        }

        new Timer().schedule(task as TimerTask, duration.toMillis())
        return target
    }
}
