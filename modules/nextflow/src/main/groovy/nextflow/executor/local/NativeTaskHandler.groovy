/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.executor.local

import java.lang.reflect.InvocationTargetException
import java.util.concurrent.Callable
import java.util.concurrent.Future

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.executor.Executor
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus

/**
 * Executes a native piece of groovy code
 */
@Slf4j
@CompileStatic
class NativeTaskHandler extends TaskHandler {

    Future<Object> result

    private Session session

    private Executor executor

    private class TaskSubmit implements Callable {

        final TaskRun task

        TaskSubmit( TaskRun obj ) { task = obj }

        @Override
        Object call() throws Exception {
            try  {
                return task.code.call()
            }
            catch( Throwable error ) {
                return error
            }
            finally {
                executor.getTaskMonitor().signal()
            }
        }
    }

    protected NativeTaskHandler(TaskRun task, Executor executor) {
        super(task)
        this.executor = executor
        this.session = executor.session
    }


    @Override
    void submit() {
        // submit for execution by using session executor service
        // it returns an error when everything is OK
        // of the exception throw in case of error
        result = session.getExecService().submit(new TaskSubmit(task))
        status = TaskStatus.SUBMITTED
    }

    @Override
    boolean checkIfRunning() {
        if( isSubmitted() && result != null ) {
            status = TaskStatus.RUNNING
            return true
        }

        return false
    }

    @Override
    boolean checkIfCompleted() {
        if( isRunning() && result.isDone() ) {
            status = TaskStatus.COMPLETED
            final ret = result.get()
            if( ret instanceof InvocationTargetException ) {
                task.error = ret.cause
            }
            else if( ret instanceof Throwable ) {
                task.error = (Throwable)ret
            }
            else {
                task.stdout = ret
            }
            return true
        }
        return false
    }

    @Override
    protected void killTask() {
        if( result ) result.cancel(true)
    }

}
