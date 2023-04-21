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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.ArrayTaskPollingMonitor
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
/**
 * Executor that submits tasks in batches to a target executor
 * that supports array jobs.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ArrayExecutor extends Executor {

    private Executor target

    @Override
    TaskMonitor createTaskMonitor() {
        final targetName = session.getExecConfigProp('array', 'target', 'local') as String
        target = session.executorFactory.getExecutor(targetName, session)

        if( target !instanceof ArrayTaskAware )
            throw new IllegalArgumentException("Executor '${targetName}' does not support array jobs")

        return ArrayTaskPollingMonitor.create(session, target, 100, Duration.of('5 sec'), 100)
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        target.createTaskHandler(task)
    }

}
