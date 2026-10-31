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
package groovyx.gpars.dataflow.operator;

import groovyx.gpars.group.PGroup;

/**
 * Runs a dataflow processor: reads its inputs, fires its body and handles its life cycle
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
abstract class OperatorEngine {

    /**
     * Starts the operator, on the given group when not null
     */
    abstract void start(PGroup group);

    abstract void terminate();

    abstract void terminateAfterNextRun();

    abstract void join() throws InterruptedException;

    /**
     * Whether the operator was started and has not stopped yet
     */
    abstract boolean isActive();

    abstract void resume();

    abstract boolean isPaused();

    abstract long getPauseCount();
}
