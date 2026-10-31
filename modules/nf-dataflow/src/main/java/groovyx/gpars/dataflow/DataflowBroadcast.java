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
package groovyx.gpars.dataflow;

import groovyx.gpars.dataflow.stream.DataflowStreamReadAdapter;
import groovyx.gpars.dataflow.stream.DataflowStreamWriteAdapter;

/**
 * A broadcast channel: each read channel receives every value bound after its creation.
 * <p>
 * Rewrites the GPars 1.2.1 class of the same name (Copyright 2008-2011 the original author or authors).
 *
 * @param <T> The type of messages to pass through the stream
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public final class DataflowBroadcast<T> extends DataflowStreamWriteAdapter<T> {

    public DataflowBroadcast() {
    }

    /**
     * Creates a read channel for the values bound from now on
     */
    public DataflowReadChannel<T> createReadChannel() {
        return new DataflowStreamReadAdapter<T>(this);
    }
}
