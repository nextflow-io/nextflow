// GPars - Groovy Parallel Systems
//
// Copyright © 2008-2012  The original author or authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Modified by Seqera Labs for Nextflow (2026): removed the select and selector factory methods, tasks, splitters and whenAllBound; the default group runs on the dataflow scheduler pool.

package groovyx.gpars.dataflow;

import groovy.lang.Closure;
import groovy.lang.DelegatesTo;
import groovyx.gpars.dataflow.operator.DataflowOperator;
import groovyx.gpars.dataflow.operator.DataflowProcessor;
import groovyx.gpars.group.DefaultPGroup;
import groovyx.gpars.group.PGroup;
import groovyx.gpars.scheduler.Scheduler;

import java.util.List;
import java.util.Map;

/**
 * Contains factory methods to create dataflow operators and start them.
 *
 * @author Vaclav Pech, Dierk Koenig
 *         Date: Jun 4, 2009
 */
@SuppressWarnings({"rawtypes", "RawUseOfParameterizedType", "AbstractClassWithoutAbstractMethods", "AbstractClassNeverImplemented", "ConstantDeclaredInAbstractClass", "UtilityClass", "unchecked"})
public abstract class Dataflow {

    /**
     * The parallel group used by all dataflow operators by default.
     */
    public static final PGroup DATA_FLOW_GROUP = new DefaultPGroup(Scheduler.pool());

    /**
     * Maps threads/tasks to parallel groups they belong to
     */
    public static final ThreadLocal<PGroup> activeParallelGroup = new ThreadLocal<PGroup>();

    /**
     * Retrieves the thread-local value of the active PGroup or the default DataflowGroup
     *
     * @return The PGroup to use for DF within the current thread
     */
    public static PGroup retrieveCurrentDFPGroup() {
        PGroup pGroup = activeParallelGroup.get();
        if (pGroup == null) {
            pGroup = Dataflow.DATA_FLOW_GROUP;
        }
        return pGroup;
    }

    /**
     * Sets the supplied PGroup as the default for the given block of code. All dataflow functions, such as task or operator
     * and callback handlers, will use the PGroup and its thread pool for their scheduling.
     *
     * @param group The group to make the default inside the block
     * @param code  The code to run with overriden default
     * @return The value returned from the supplied code block
     */
    public static Object usingGroup(final PGroup group, final Closure code) {
        final PGroup original = activeParallelGroup.get();
        try {
            activeParallelGroup.set(group);
            return code.call();
        } finally {
            activeParallelGroup.set(original);
        }
    }

    /**
     * Creates an operator using the default dataflow parallel group
     *
     * @param channels A map specifying "inputs" and "outputs" - dataflow channels (instances of the DataflowQueue or DataflowVariable classes) to use for inputs and outputs
     * @param code     The operator's body to run each time all inputs have a value to read
     * @return A new active operator instance
     */
    public static DataflowProcessor operator(final Map channels, @DelegatesTo(DataflowOperator.class) final Closure code) {
        return retrieveCurrentDFPGroup().operator(channels, code);
    }

    /**
     * Creates an operator using the current parallel group
     *
     * @param inputChannels  dataflow channels to use for input
     * @param outputChannels dataflow channels to use for output
     * @param code           The operator's body to run each time all inputs have a value to read
     * @return A new active operator instance
     */
    public static DataflowProcessor operator(final List inputChannels, final List outputChannels, @DelegatesTo(DataflowOperator.class) final Closure code) {
        return retrieveCurrentDFPGroup().operator(inputChannels, outputChannels, code);
    }

    /**
     * Creates an operator using the current parallel group
     *
     * @param inputChannels  dataflow channels to use for input
     * @param outputChannels dataflow channels to use for output
     * @param maxForks       Number of parallel threads running operator's body, defaults to 1
     * @param code           The operator's body to run each time all inputs have a value to read
     * @return A new active operator instance
     */
    public static DataflowProcessor operator(final List inputChannels, final List outputChannels, final int maxForks, @DelegatesTo(DataflowOperator.class) final Closure code) {
        return retrieveCurrentDFPGroup().operator(inputChannels, outputChannels, maxForks, code);
    }

    /**
     * Creates an operator using the current parallel group
     *
     * @param input  a dataflow channel to use for input
     * @param output a dataflow channel to use for output
     * @param code   The operator's body to run each time all inputs have a value to read
     * @return A new active operator instance
     */
    public static DataflowProcessor operator(final DataflowReadChannel input, final DataflowWriteChannel output, @DelegatesTo(DataflowOperator.class) final Closure code) {
        return retrieveCurrentDFPGroup().operator(input, output, code);
    }

    /**
     * Creates an operator using the current parallel group
     *
     * @param input    a dataflow channel to use for input
     * @param output   a dataflow channel to use for output
     * @param maxForks Number of parallel threads running operator's body, defaults to 1
     * @param code     The operator's body to run each time all inputs have a value to read
     * @return A new active operator instance
     */
    public static DataflowProcessor operator(final DataflowReadChannel input, final DataflowWriteChannel output, final int maxForks, @DelegatesTo(DataflowOperator.class) final Closure code) {
        return retrieveCurrentDFPGroup().operator(input, output, maxForks, code);
    }

}
