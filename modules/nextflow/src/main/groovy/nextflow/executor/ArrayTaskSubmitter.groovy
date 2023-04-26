/*
 * Copyright 2013-2023, Seqera Labs
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

import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.TaskHandler

/**
 * Submit tasks as an array job.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ArrayTaskSubmitter {

    protected List<TaskHandler> array

    private AtomicInteger collected = new AtomicInteger()

    ArrayTaskSubmitter(List<TaskHandler> array) {
        this.array = array

        for( TaskHandler handler : array )
            handler.arraySubmitter = this
    }

    List<TaskHandler> getArray() { array }

    /**
     * Mark a task as ready to be submitted.
     *
     * @param handler
     */
    void collect(TaskHandler handler) {
        if( collected.incrementAndGet() == array.size() )
            submit()
    }

    /**
     * Submit the array job.
     *
     * By default, this method simply submits each task individually.
     * It should be overridden to submit an array job to the underlying
     * executor.
     */
    protected void submit() {
        for( TaskHandler handler : array )
            handler.submit()
    }

}
