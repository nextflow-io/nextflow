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

package nextflow.dataflow;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import groovyx.gpars.scheduler.Scheduler;

/**
 * Runs calls on the dataflow pool, e.g. on the loop thread with {@code NXF_DATAFLOW_POOL_SIZE=1}.
 */
public final class DataflowThreads {

    private DataflowThreads() {}

    /**
     * Submits the call {@code n} times to the dataflow pool and returns the results in order
     *
     * @throws Exception the error of the first failed call
     */
    public static <T> List<T> call(int n, Callable<T> call) throws Exception {
        final List<Future<T>> futures = new ArrayList<>(n);
        for( int i = 0; i < n; i++ )
            futures.add(Scheduler.pool().submit(call));
        final List<T> result = new ArrayList<>(n);
        for( Future<T> future : futures ) {
            try {
                result.add(future.get(30, TimeUnit.SECONDS));
            }
            catch( ExecutionException e ) {
                if( e.getCause() instanceof Exception cause )
                    throw cause;
                throw e;
            }
        }
        return result;
    }
}
