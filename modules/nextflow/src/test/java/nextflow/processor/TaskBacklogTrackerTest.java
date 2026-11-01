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

package nextflow.processor;

import java.util.List;

import groovyx.gpars.dataflow.DataflowQueue;
import groovyx.gpars.dataflow.operator.DataflowProcessor;
import nextflow.dataflow.engine.Backpressure;
import nextflow.dataflow.engine.ProcessBacklog;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TaskBacklogTrackerTest {

    private ProcessBacklog backlog;

    private TaskBacklogTracker tracker;

    @BeforeEach
    void setup() {
        // the last input is the control channel
        final var inputs = List.of(new DataflowQueue<>(), new DataflowQueue<>());
        final var operator = new DataflowProcessor(inputs, null, List.of(), 1, new InvokeTaskAdapter(null, 2));
        backlog = new Backpressure(null, 4).register("foo", inputs, operator);
        tracker = new TaskBacklogTracker(() -> backlog);
    }

    @Test
    void shouldMoveAScheduledTaskFromPreparingToPending() {
        final var task = new TaskRun();

        tracker.accepted(task);
        assertThat(backlog.getPreparing()).isEqualTo(1);

        tracker.scheduled(task);
        assertThat(backlog.getPreparing()).isZero();
        assertThat(backlog.getPending()).isEqualTo(1);

        // the end of the preparation does not release a scheduled task
        tracker.prepared(task);
        assertThat(backlog.getPreparing()).isZero();
        assertThat(backlog.getPending()).isEqualTo(1);

        tracker.dequeued(task);
        assertThat(backlog.getPending()).isZero();
    }

    @Test
    void shouldReleaseADiscardedTask() {
        final var task = new TaskRun();

        tracker.accepted(task);
        tracker.prepared(task);
        tracker.prepared(task);

        assertThat(backlog.getPreparing()).isZero();
        assertThat(backlog.getPending()).isZero();
        assertThat(tracker.getPreparing()).isZero();
    }

    @Test
    void shouldCountAnAcceptedTaskOnce() {
        final var task = new TaskRun();

        tracker.accepted(task);
        tracker.accepted(task);
        assertThat(backlog.getPreparing()).isEqualTo(1);

        tracker.prepared(task);
        assertThat(backlog.getPreparing()).isZero();
    }

    @Test
    void shouldCountATaskScheduledWithoutPreparationAsPending() {
        // e.g. a retried task or a task array
        final var task = new TaskRun();

        tracker.scheduled(task);
        tracker.scheduled(task);
        assertThat(backlog.getPreparing()).isZero();
        assertThat(backlog.getPending()).isEqualTo(1);

        tracker.dequeued(task);
        tracker.dequeued(task);
        assertThat(backlog.getPending()).isZero();
        assertThat(tracker.getPending()).isZero();
    }

    @Test
    void shouldIgnoreATaskThatWasNotCounted() {
        tracker.prepared(new TaskRun());
        tracker.dequeued(new TaskRun());

        assertThat(backlog.getPreparing()).isZero();
        assertThat(backlog.getPending()).isZero();
    }

    @Test
    void shouldCountAnArrayChildAsPreparedAndTheArrayAsPending() {
        final var child = new TaskRun();
        final var array = new TaskRun();

        tracker.accepted(child);
        // the array is scheduled while its last child is still being prepared
        tracker.scheduled(array);
        tracker.prepared(child);

        assertThat(backlog.getPreparing()).isZero();
        assertThat(backlog.getPending()).isEqualTo(1);

        tracker.dequeued(array);
        assertThat(backlog.getPending()).isZero();
    }

    @Test
    void shouldNotTrackTasksWithoutBackpressure() {
        final var tracker = new TaskBacklogTracker(() -> null);
        final var task = new TaskRun();

        tracker.accepted(task);
        tracker.scheduled(task);

        assertThat(tracker.getPreparing()).isZero();
        assertThat(tracker.getPending()).isZero();
    }
}
