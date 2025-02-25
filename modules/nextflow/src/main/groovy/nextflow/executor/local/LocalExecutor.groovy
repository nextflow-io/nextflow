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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.Executor
import nextflow.executor.SupportedScriptTypes
import nextflow.fusion.FusionHelper
import nextflow.processor.LocalPollingMonitor
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
/**
 * Executes the specified task on the locally exploiting the underlying Java thread pool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@SupportedScriptTypes( [ScriptType.SCRIPTLET, ScriptType.GROOVY] )
class LocalExecutor extends Executor {

    @Override
    protected TaskMonitor createTaskMonitor() {
        return LocalPollingMonitor.create(session, name)
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir

        if( task.type == ScriptType.GROOVY )
            return new NativeTaskHandler(task,this)
        else
            return new LocalTaskHandler(task,this)
    }

    @Override
    boolean isContainerNative() {
        return isFusionEnabled()
    }

    @Override
    boolean isFusionEnabled() {
        return FusionHelper.isFusionEnabled(session)
    }
}

