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
package nextflow.dataflow.engine;

import java.util.List;

import groovyx.gpars.dataflow.DataflowBroadcast;
import groovyx.gpars.dataflow.DataflowQueue;
import groovyx.gpars.dataflow.DataflowReadChannel;
import groovyx.gpars.dataflow.DataflowVariable;
import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ProcessBacklogTest {

    @SuppressWarnings("rawtypes")
    private static ProcessBacklog backlog(Backpressure bp, DataflowReadChannel... inputs) {
        return new ProcessBacklog(bp, "foo", List.of(inputs), null);
    }

    private static void times(int n, Runnable action) {
        for( int i = 0; i < n; i++ )
            action.run();
    }

    @Test
    void shouldCloseTheOperatorGateWhilePreparingAndPendingTasksReachTheHighWaterMark() {
        var backlog = backlog(new Backpressure(null, 8));

        times(3, backlog::taskStarted);
        times(5, backlog::taskPending);
        assertThat(backlog.getPreparing()).isEqualTo(3);
        assertThat(backlog.getPending()).isEqualTo(5);
        assertThat(backlog.getAsBoolean()).isTrue();

        backlog.taskPrepared();
        assertThat(backlog.getAsBoolean()).isFalse();

        backlog.taskPending();
        backlog.close();
        assertThat(backlog.getAsBoolean()).isFalse();
    }

    @Test
    void shouldOpenTheOperatorGateWhenBackpressureIsReleased() {
        var bp = new Backpressure(null, 2);
        var backlog = backlog(bp);
        times(2, backlog::taskPending);
        assertThat(backlog.getAsBoolean()).isTrue();

        bp.release();
        assertThat(backlog.getAsBoolean()).isFalse();
    }

    @Test
    void shouldCountTheCompleteInputSets() {
        var a = new DataflowBroadcast<Object>();
        var b = new DataflowBroadcast<Object>();
        var q = new DataflowQueue<Object>();
        var backlog = backlog(new Backpressure(null, 8), a.createReadChannel(), b.createReadChannel(), q);
        assertThat(backlog.ready(100)).isZero();

        for( int i = 0; i < 5; i++ ) a.bind(i);
        for( int i = 0; i < 2; i++ ) b.bind(i);
        for( int i = 0; i < 3; i++ ) q.bind(i);
        assertThat(backlog.ready(100)).isEqualTo(2);
        assertThat(backlog.ready(1)).isEqualTo(1);
    }

    @Test
    void shouldNotCountItemsWaitingForAnUnboundValue() {
        var ch = new DataflowBroadcast<Object>();
        var value = new DataflowVariable<Object>();
        var backlog = backlog(new Backpressure(null, 4), ch.createReadChannel(), value);

        for( int i = 0; i < 10; i++ ) ch.bind(i);
        assertThat(backlog.ready(100)).isZero();
        assertThat(backlog.isOver(4)).isFalse();

        value.bind("x");
        assertThat(backlog.ready(100)).isEqualTo(10);
        assertThat(backlog.isOver(4)).isTrue();
    }

    @Test
    void shouldNotCountAProcessWithOnlyValueInputs() {
        var value = new DataflowVariable<Object>();
        value.bind(1);
        var backlog = backlog(new Backpressure(null, 4), value);

        assertThat(backlog.ready(100)).isZero();
    }

    @Test
    void shouldCompareTheWholeBacklogWithAMark() {
        var ch = new DataflowBroadcast<Object>();
        var backlog = backlog(new Backpressure(null, 8), ch.createReadChannel());

        times(2, backlog::taskStarted);
        times(3, backlog::taskPending);
        ch.bind(0);
        ch.bind(1);
        assertThat(backlog.isOver(8)).isFalse();
        assertThat(backlog.isOver(7)).isTrue();

        backlog.close();
        assertThat(backlog.isOver(1)).isFalse();
    }
}
