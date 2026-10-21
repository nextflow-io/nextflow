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

package nextflow.executor

import java.util.concurrent.Callable
import java.util.concurrent.ExecutionException
import java.util.concurrent.ExecutorService
import java.util.concurrent.Future
import java.util.concurrent.FutureTask
import java.util.concurrent.LinkedBlockingQueue

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.TaskHandler
import nextflow.util.ThreadPoolBuilder
/**
 * Bounded thread pool that prepares the launcher scripts of job array
 * children in parallel, for executors whose work directory is on remote
 * storage (where each launcher requires several remote writes).
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ArrayLauncherPool {

    static final int DEFAULT_SIZE = 16

    private final String name

    private final int size

    private ExecutorService pool

    ArrayLauncherPool(String name, int size=DEFAULT_SIZE) {
        this.name = name
        this.size = size
    }

    /**
     * Prepare the launcher of each child task, failing as soon as
     * any child fails.
     *
     * @param children
     */
    void prepare(List<TaskHandler> children) {
        final executor = getPool()
        final completed = new LinkedBlockingQueue<Future<Object>>()
        final futures = new ArrayList<Future<Object>>(children.size())
        try {
            for( final child : children ) {
                final task = new NotifyingTask(prepare0(child), completed)
                futures.add(task)
                executor.execute(task)
            }
            for( int i=0; i<futures.size(); i++ )
                completed.take().get()
        }
        catch( ExecutionException e ) {
            throw e.cause
        }
        catch( InterruptedException e ) {
            Thread.currentThread().interrupt()
            throw e
        }
        finally {
            for( final future : futures )
                future.cancel(true)
        }
    }

    private Callable<Object> prepare0(TaskHandler child) {
        return {
            try {
                child.prepareLauncher()
                return null
            }
            catch( Throwable e ) {
                log.warn "Unable to prepare launcher for job array task `${child.task.lazyName()}` -- work-dir=${child.task.getWorkDirStr()}"
                throw e
            }
        } as Callable<Object>
    }

    /**
     * Future that notifies the caller when it is done, including when
     * it is cancelled while still queued.
     */
    private static class NotifyingTask extends FutureTask<Object> {

        private final Queue<Future<Object>> completed

        NotifyingTask(Callable<Object> callable, Queue<Future<Object>> completed) {
            super(callable)
            this.completed = completed
        }

        @Override
        protected void done() {
            completed.add(this)
        }
    }

    synchronized protected ExecutorService getPool() {
        if( pool == null ) {
            pool = new ThreadPoolBuilder()
                .withName(name)
                .withMinSize(size)
                .withMaxSize(size)
                .withAllowCoreThreadTimeout(true)
                .build()
        }
        return pool
    }

    /**
     * Shutdown the pool, cancelling any pending preparations so that
     * no caller remains blocked waiting on them.
     */
    synchronized void shutdown() {
        if( pool == null )
            return
        for( final task : pool.shutdownNow() ) {
            if( task instanceof Future )
                ((Future) task).cancel(true)
        }
    }

}
