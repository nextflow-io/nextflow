/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.ChoiceClosure
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.Global
import nextflow.Session
/**
 * Implements the logic for {@link DataflowExtensions#choice} operator(s)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ChoiceOp {

    /**
     * The operator source channel "normalised" in a list object.
     * It must contain exactly *one* DataflowReadChannel instance
     */
    private List<DataflowReadChannel> source

    /**
     * The list of output channels resulting from the choice operation
     */
    private List<DataflowWriteChannel> outputs

    /**
     * A closure implementing the *choice* strategy. It returns the index of the
     * selected channel given the actual item emitted by the source channel
     */
    private Closure<Integer> code

    /**
     * {@code true} when the `choice` is applied to a {@link groovyx.gpars.dataflow.DataflowVariable}
     */
    private boolean stopOnFirst

    /**
     * The current nextflow {@link Session} object
     */
    private Session session = (Session)Global.session

    /**
     * Creates the choice operator
     *
     * @param source The source channel either a {@link groovyx.gpars.dataflow.DataflowQueue} or a {@link groovyx.gpars.dataflow.DataflowVariable}
     * @param outputs The list of output channels
     * @param code The closure implementing the *choice* strategy. See {@link #code}
     */
    ChoiceOp(DataflowReadChannel source, List<DataflowWriteChannel> outputs, Closure<Integer> code) {
        assert source
        assert outputs
        this.source = [source]
        this.outputs = outputs
        this.code = code
        this.stopOnFirst = source instanceof DataflowExpression
    }

    /**
     * @return A {@link DataflowEventAdapter} that close properly the output
     * channels when required
     */
    private createListener() {

        def result = new DataflowEventAdapter() {
            @Override
            void afterRun(DataflowProcessor processor, List<Object> messages) {
                if( !stopOnFirst ) return
                // -- terminate the process
                processor.terminate()
                // -- close the output channels
                outputs.each {
                    if( !(it instanceof DataflowExpression))
                        it.bind(Channel.STOP)

                    else if( !(it as DataflowExpression).isBound() )
                        it.bind(Channel.STOP)
                }
            }

            @Override
            public boolean onException(final DataflowProcessor processor, final Throwable e) {
                log.error("@unknown", e)
                session.abort(e)
                return true;
            }
        }

        return [result]
    }

    /**
     * Applies the choice operator
     */
    def apply() {

        def params = [
                inputs: source,
                outputs: outputs,
                listeners: createListener()
        ]

        DataflowHelper.newOperator(params, new ChoiceClosure(code))
    }


}
