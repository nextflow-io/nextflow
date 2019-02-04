/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.processor

/**
 * Declares the methods for scheduling and handling tasks
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface TaskMonitor {

    /**
     * Schedule a {@code TaskHandler} instance adding it to the queue of tasks to be processed.
     * Subclasses provide concrete logic to manage them.
     *
     * Note: depending the task scheduling implementation this operation may be blocking
     * to await there's free space in the tasks queue
     *
     * @param handler {@code TaskHandler} instance
     */
    void schedule(TaskHandler handler)

    /**
     * Remove the {@code TaskHandler} instance from the queue of tasks to be processed
     *
     * @param handler A not null {@code TaskHandler} instance
     */
    boolean evict(TaskHandler handler)

    /**
     * Start the monitoring activity for the queued tasks
     * @return The instance itself, useful to chain methods invocation
     */
    TaskMonitor start()

    /**
     * Notify when a task terminates
     */
    void signal()
}
