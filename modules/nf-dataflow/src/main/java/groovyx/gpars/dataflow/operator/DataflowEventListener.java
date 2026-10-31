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
// Modified by Seqera Labs for Nextflow (2026): merged DataflowEventAdapter as default methods; removed the
// registered, afterStart, messageSentOut and customEvent events.

package groovyx.gpars.dataflow.operator;

import groovyx.gpars.dataflow.DataflowReadChannel;

import java.util.List;

/**
 * Enables external entities to observe the state of dataflow operators. Every event has a default no-op handler.
 *
 * @author Vaclav Pech
 */
public interface DataflowEventListener {

    /**
     * Invoked immediately after the operator terminates
     *
     * @param processor The reporting dataflow operator
     */
    default void afterStop(DataflowProcessor processor) {
    }

    /**
     * Invoked if an exception occurs.
     * If any of the listeners returns true, the operator will terminate.
     * Exceptions outside of the operator's body will terminate the operator irrespective of the listeners' votes.
     * When using maxForks, the method may be invoked from threads running the forks.
     *
     * @param processor The reporting dataflow operator
     * @param e         The thrown exception
     * @return True, if the operator should terminate in response to the exception, false otherwise. True by default.
     */
    default boolean onException(DataflowProcessor processor, Throwable e) {
        return true;
    }

    /**
     * Invoked when a message becomes available in an input channel.
     *
     * @param processor The reporting dataflow operator
     * @param channel   The input channel holding the message
     * @param index     The index of the input channel within the operator
     * @param message   The incoming message
     * @return The original message or a message that should be used instead
     */
    default Object messageArrived(DataflowProcessor processor, DataflowReadChannel<Object> channel, int index, Object message) {
        return message;
    }

    /**
     * Invoked when a control message (a {@link PoisonPill}) becomes available in an input channel.
     *
     * @param processor The reporting dataflow operator
     * @param channel   The input channel holding the message
     * @param index     The index of the input channel within the operator
     * @param message   The incoming message
     * @return The original message or a message that should be used instead
     */
    default Object controlMessageArrived(DataflowProcessor processor, DataflowReadChannel<Object> channel, int index, Object message) {
        return message;
    }

    /**
     * Invoked when all messages required to trigger the operator become available in the input channels.
     *
     * @param processor The reporting dataflow operator
     * @param messages  The incoming messages
     * @return The original list of messages or a modified/new list of messages that should be used instead
     */
    default List<Object> beforeRun(DataflowProcessor processor, List<Object> messages) {
        return messages;
    }

    /**
     * Invoked when the operator completes a single run
     * When using maxForks, the method may be invoked from threads running the forks.
     *
     * @param processor The reporting dataflow operator
     * @param messages  The incoming messages that have been processed
     */
    default void afterRun(DataflowProcessor processor, List<Object> messages) {
    }
}
