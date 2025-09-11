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

import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowEventListener
import groovyx.gpars.dataflow.operator.DataflowOperator
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.Global
import nextflow.Session
import nextflow.dag.NodeMarker
import nextflow.extension.DataflowHelper
import nextflow.prov.OperatorRun
import nextflow.prov.Prov
import nextflow.prov.Tracker
/**
 * Model a dataflow operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class Op {

    /**
     * Associate a {@link DataflowProcessor} with the corresponding {@link OpContext} object
     * tracking running operator context
     */
    static final public ConcurrentHashMap<DataflowProcessor, OpContext> allContexts = new ConcurrentHashMap<>();

    static List<Object> unwrap(List messages) {
        final ArrayList<Object> result = new ArrayList<>();
        for( Object it : messages ) {
            result.add(it instanceof Tracker.Msg ? ((Tracker.Msg)it).value : it);
        }
        return result;
    }

    static Object unwrap(Object it) {
        return it instanceof Tracker.Msg ? it.value : it
    }

    static Tracker.Msg wrap(Object obj) {
        obj instanceof Tracker.Msg ? obj : Tracker.Msg.of(obj)
    }

    static void bind(DataflowProcessor dp, DataflowWriteChannel channel, Object msg) {
        if( dp==null )
            throw new IllegalStateException("DataflowProcessor argument cannot be null")
        try {
            if( msg instanceof PoisonPill ) {
                final emit = channel !instanceof DataflowExpression || !(channel as DataflowChannel).isBound()
                if( emit )
                    channel.bind(msg)
                allContexts.remove(dp)
            }
            else {
                final ctx = allContexts.get(dp)
                if( !ctx )
                    throw new IllegalStateException("Cannot find any context for operator=$dp")
                final run = ctx.getOperatorRun()
                Prov.getTracker().bindOutput(run, channel, msg)
            }
        }
        catch (Throwable t) {
            log.error("Unexpected resolving execution provenance: ${t.message}", t)
            (Global.session as Session).abort(t)
        }
    }

    /**
     * This method gets a list of one or more {@link OpDatum} objects, and binds
     * the corresponding values to the specified channel using a new {@link OperatorRun}
     * instance using as {@link OperatorRun#inputIds} the set of IDs obtained by the
     * provide "entries" lists.
     *
     * @param target
     *      The dataflow channel where outputs are going to be bound.
     * @param entries
     *      A list of one or more {@link OpDatum} holding the value to be bound.
     * @param singleton
     *      When {@code true} and the entries list is provided, the value object is emitted
     *      in place of a list with a sole object, when {@code false} the list of values
     *      is emitted in any case.
     */
    static void bindRunValues(DataflowWriteChannel target, List<OpDatum> entries, boolean singleton) {
        final inputs = new ArrayList<Integer>(entries.size())
        final values = new ArrayList(entries.size())
        for( Object it : entries ) {
            if( it instanceof OpDatum ) {
                inputs.addAll(it.run.inputIds)
                values.add(it.value)
            }
            else
                values.add(it)
        }
        final run = new OperatorRun(new LinkedHashSet<Integer>(inputs))
        final out = singleton && values.size()==1 ? values[0] : values
        Op.bind(run, target, out)
    }

    static void bind(OperatorRun run, DataflowWriteChannel channel, Object msg) {
        Prov.getTracker().bindOutput(run, channel, msg)
    }

    private static Session getSession() { Global.getSession() as Session }

    private List<DataflowReadChannel> inputs
    private List<DataflowWriteChannel> outputs
    private List<DataflowEventListener> listeners
    private OpContext context = new ContextSequential()
    private Closure code

    List<DataflowReadChannel> getInputs() { inputs }
    List<DataflowWriteChannel> getOutputs() { outputs }
    List<DataflowEventListener> getListeners() { listeners }
    OpContext getContext() { context }
    Closure getCode() { code }

    /**
     * Set the operator input channel
     *
     * @param channel
     *      The channel from which the operator is receiving the input data.
     * @return
     *      The operator object itself.
     */
    Op withInput(DataflowReadChannel channel) {
        assert channel != null
        this.inputs = List.of(channel)
        return this
    }

    /**
     * Set the operator input channels
     *
     * @param channels
     *      The channels from which the operator is receiving the input data
     * @return
     *      The operator object itself.
     */
    Op withInputs(List<DataflowReadChannel> channels) {
        assert channels != null
        this.inputs = channels
        return this
    }

    /**
     * Set the operator output channel
     *
     * @param channel
     *      The channel to which the operator is binding output data.
     * @return
     *      The operator object itself.
     */
    Op withOutput(DataflowWriteChannel channel) {
        assert channel != null
        this.outputs = List.of(channel)
        return this
    }

    /**
     * Set the operator output channels
     *
     * @param channels
     *      The channels to which the operator is binding output data.
     * @return
     *      The operator object itself.
     */
    Op withOutputs(List<DataflowWriteChannel> channels) {
        assert channels != null
        this.outputs = channels
        return this
    }

    /**
     * Set a {@link DataflowEventListener} associated with the operator.
     *
     * @param listener
     *      The {@link DataflowEventListener} listener object.
     * @return
     *      The operator object itself.
     */
    Op withListener(DataflowEventListener listener) {
        if( listener )
            this.listeners = List.of(listener)
        return this
    }

    /**
     * Set one or more {@link DataflowEventListener} objects associated with the operator.
     *
     * @param listeners
     *      The list of {@link DataflowEventListener} listener objects.
     * @return
     *      The operator object itself.
     */
    Op withListeners(List<DataflowEventListener> listeners) {
        if( listeners )
            this.listeners = listeners
        return this
    }

    /**
     * Set the operator inputs, outputs and listeners parameters using a map object.
     *
     * @param params
     *      A {@link Map} object holding the {@code inputs} and {@code outputs} channels
     *      and {@code listeners} objects.
     * @return
     *      The operator object itself.
     */
    Op withParams(Map params) {
        if( params.inputs )
            this.inputs = params.inputs as List<DataflowReadChannel>
        if( params.outputs )
            this.outputs = params.outputs as List<DataflowWriteChannel>
        if( params.listeners )
            this.listeners = params.listeners as List<DataflowEventListener>
        return this
    }

    /**
     * Set the context object associated associated with the operator run.
     *
     * @param context
     *      The {@link OpContext} object associated with the operator execution.
     * @return
     *      The operator object itself.
     */
    Op withContext(OpContext context) {
        if( context!=null )
            this.context = context
        return this
    }

    /**
     * Set the operator action as a closure
     *
     * @param code
     *      The action that should be performed when the operator run.
     * @return
     *      The operator object itself.
     *
     */
    Op withCode(Closure code) {
        this.code = code
        return this
    }

    /**
     * Map the operator inputs, outputs and listeners to a Java {@Link Map}
     * for compatibility with {@link DataflowProcessor} API
     *
     * @return
     *      A map object holding the operator {@code inputs}, {@code outputs} channels
     *      and the {@code listeners} object.
     */
    protected Map toMap() {
        final ret = new HashMap()
        ret.inputs = inputs ?: List.of()
        ret.outputs = outputs ?: List.of()
        ret.listeners = listeners ?: List.of()
        return ret
    }

    /**
     * Creates the {@link DataflowOperator} object and start it
     *
     * @return
     *      The corresponding {@link DataflowOperator} object
     */
    DataflowProcessor apply() {
        assert inputs
        assert code
        assert context

        if( !listeners )
            listeners = List.of(DataflowHelper.defaultErrorListener())

        // Encapsulate the target "code" closure with a "OpClosure" object
        // to grab input data and track the execution provenance
        final closure = new OpClosure(code, context)
        final group = Dataflow.retrieveCurrentDFPGroup()
        final operator = new DataflowOperator(group, toMap(), closure)
        allContexts.put(operator, context)
        operator.start()
        // track the operator as dag node
        NodeMarker.appendOperator(operator)
        if( session && session.allOperators != null ) {
            session.allOperators << operator
        }
        return operator
    }

}
