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

package nextflow.extension

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.Global
import nextflow.Session
import nextflow.prov.OperatorRun
import nextflow.prov.Prov
import nextflow.prov.Tracker

/**
 * Operator helpers methods
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class Op {

    static final @PackageScope ThreadLocal<OperatorRun> currentOperator = new ThreadLocal<>()

    static List<Object> unwrap(List messages) {
        return messages.collect(it -> it instanceof Tracker.Msg ? it.value : it)
    }

    static Object unwrap(Object it) {
        return it instanceof Tracker.Msg ? it.value : it
    }

    static Tracker.Msg wrap(Object obj) {
        obj instanceof Tracker.Msg ? obj : Tracker.Msg.of(obj)
    }

    static void bind(DataflowWriteChannel channel, Object msg) {
        try {
            if( msg instanceof PoisonPill )
                channel.bind(msg)
            else
                Prov.getTracker().bindOutput(currentOperator.get(), channel, msg)
        }
        catch (Throwable t) {
            log.error("Unexpected resolving execution provenance: ${t.message}", t)
            (Global.session as Session).abort(t)
        }
    }

    static Closure instrument(Closure op, boolean accumulator=false) {
        return new InvokeOperatorAdapter(op, accumulator)
    }

    static class InvokeOperatorAdapter extends Closure {

        private final Closure target

        private final boolean accumulator

        private OperatorRun previousRun

        private InvokeOperatorAdapter(Closure code, boolean accumulator) {
            super(code.owner, code.thisObject)
            this.target = code
            this.target.delegate = code.delegate
            this.target.setResolveStrategy(code.resolveStrategy)
            this.accumulator = accumulator
        }

        @Override
        Class<?>[] getParameterTypes() {
            return target.getParameterTypes()
        }

        @Override
        int getMaximumNumberOfParameters() {
            return target.getMaximumNumberOfParameters()
        }

        @Override
        Object getDelegate() {
            return target.getDelegate()
        }

        @Override
        Object getProperty(String propertyName) {
            return target.getProperty(propertyName)
        }

        @Override
        int getDirective() {
            return target.getDirective()
        }

        @Override
        void setDelegate(Object delegate) {
            target.setDelegate(delegate)
        }

        @Override
        void setDirective(int directive) {
            target.setDirective(directive)
        }

        @Override
        void setResolveStrategy(int resolveStrategy) {
            target.setResolveStrategy(resolveStrategy)
        }

        @Override
        void setProperty(String propertyName, Object newValue) {
            target.setProperty(propertyName, newValue)
        }

        @Override
        @CompileDynamic
        Object call(final Object... args) {
            // when the accumulator flag true, re-use the previous run object
            final run = !accumulator || previousRun==null
                ? new OperatorRun()
                : previousRun
            // set as the current run in the thread local
            currentOperator.set(run)
            // map the inputs
            final inputs = Prov.getTracker().receiveInputs(run, args.toList())
            final arr = inputs.toArray()
            // todo: the spread operator should be replaced with proper array
            final ret = target.call(*arr)
            // track the previous run
            if( accumulator )
                previousRun = run
            // return the operation result
            return ret
        }

        Object call(Object args) {
            // todo: this should invoke the above one
            target.call(args)
        }

        @Override
        Object call() {
            // todo: this should invoke the above one
            target.call()
        }
    }

}
