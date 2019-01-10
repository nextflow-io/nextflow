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

package nextflow.executor

import java.nio.file.Paths

import groovy.util.logging.Slf4j
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.util.Duration

/**
 * Dummy executor, only for test purpose
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class NopeExecutor extends Executor {

    @Override
    protected TaskMonitor createTaskMonitor() {
        return TaskPollingMonitor.create(session, name, 5, Duration.of('50ms'))
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        return new NopeTaskHandler(task)
    }

}


@Slf4j
class NopeTaskHandler extends TaskHandler {

    protected NopeTaskHandler(TaskRun task) {
        super(task)
    }

    @Override
    void submit() {
        log.info ">> launching nope process: ${task}"
        task.workDir = Paths.get('.').complete()
        status = TaskStatus.SUBMITTED
        task.stdout = task.script
        task.exitStatus = 0
    }

    @Override
    boolean checkIfRunning() {
        log.debug "isRunning: $status"
        if( isSubmitted() ) {
            status = TaskStatus.RUNNING
            return true
        }
        return false
    }

    @Override
    boolean checkIfCompleted() {
        log.debug "isTerminated: $status"
        if( isRunning() ) {
            status = TaskStatus.COMPLETED
            return true
        }
        false
    }

    @Override
    void kill() { }

}

