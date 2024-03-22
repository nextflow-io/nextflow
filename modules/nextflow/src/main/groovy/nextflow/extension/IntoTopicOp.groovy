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
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.CopyChannelsClosure
import static nextflow.extension.DataflowHelper.newOperator
/**
 * Implements the {@link OperatorImpl#topic} operator
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class IntoTopicOp {

    private DataflowReadChannel source

    private String name

    private List<DataflowWriteChannel> outputs

    IntoTopicOp( DataflowReadChannel source, String name ) {
        this.source = source
        this.name = name
    }

    DataflowWriteChannel apply() {
        final target = CH.createBy(source)
        final topicSource = CH.createTopicSource(name)
        this.outputs = [target, topicSource]
        newOperator([source], outputs, new ChainWithClosure(new CopyChannelsClosure()))
        return target
    }

    List<DataflowWriteChannel> getOutputs() {
        return outputs
    }

}
