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
 * The legacy engine: runs an operator on a GPars actor ({@code NXF_DATAFLOW_ENGINE=legacy})
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
final class ActorEngine extends OperatorEngine {

    private final DataflowProcessorActor actor;

    ActorEngine(final DataflowProcessorActor actor) {
        this.actor = actor;
    }

    @Override
    void start(final PGroup group) {
        if (group != null)
            actor.setParallelGroup(group);
        actor.start();
    }

    @Override
    void terminate() {
        actor.terminate();
    }

    @Override
    void terminateAfterNextRun() {
        actor.send(StopGently.getInstance());
    }

    @Override
    void join() throws InterruptedException {
        actor.join();
    }

    @Override
    boolean isActive() {
        return actor.isActive();
    }

    @Override
    void resume() {
        actor.resume();
    }

    @Override
    boolean isPaused() {
        return actor.isPaused();
    }

    @Override
    long getPauseCount() {
        return actor.getPauseCount();
    }
}
