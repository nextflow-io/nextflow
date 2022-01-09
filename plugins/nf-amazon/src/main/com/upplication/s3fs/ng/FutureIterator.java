/*
 * Copyright 2020-2022, Seqera Labs
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

package com.upplication.s3fs.ng;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.function.Function;

/**
 * Implements an iterator that progressively submit a collection of tasks to the
 * specifies executor and iterates over the responses returned as {@link Future}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Jordi Deu-Pons <jordi@seqera.io>
 */
public class FutureIterator<REQ,RESP> implements Iterator<Future<RESP>> {

    final private ExecutorService executor;
    final private Iterator<REQ> parts;
    final private Queue<Future<RESP>> futures = new LinkedList<>();
    final private Function<REQ, RESP> task;
    final private int initialSize;

    FutureIterator(List<REQ> parts, Function<REQ, RESP> task, ExecutorService executor, int initialSize) {
        this.parts = parts.iterator();
        this.task = task;
        this.executor = executor;
        this.initialSize = initialSize;

        init();
    }

    private void init() {
        // Add up to `numWorkers` *2 parts on start
        int submitted = 0;
        while (parts.hasNext() && submitted++ < initialSize ) {
            // note: making `parts.next()` inline in the lambda causes to delay
            // the evaluate in a separate thread causing concurrency problems
            REQ req = parts.next();
            futures.add(executor.submit( () -> task.apply(req) ));
        }
    }

    @Override
    public boolean hasNext() {
        return !futures.isEmpty() || parts.hasNext();
    }

    @Override
    public Future<RESP> next() {
        // keep busy the download workers adding a new chunk
        // to download each time one is consumed
        if( parts.hasNext() ) {
            // note: making `parts.next()` inline in the lambda causes to delay
            // the evaluate in a separate thread causing concurrency problems
            REQ req = parts.next();
            futures.add(executor.submit( () -> task.apply(req)) );
        }
        try {
            return futures.poll();
        }
        catch (Throwable t) {
            // in case of error cancel all pending tasks  
            for( Future<RESP> it : futures ) {
                it.cancel(true);
            }
            throw t;
        }
    }
}
