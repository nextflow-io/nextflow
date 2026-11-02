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

package nextflow.executor.local

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.executor.Executor
import nextflow.executor.SupportedScriptTypes
import nextflow.fusion.FusionHelper
import nextflow.processor.LocalPollingMonitor
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
import nextflow.util.SysHelper
/**
 * Executes the specified task on the locally exploiting the underlying Java thread pool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@SupportedScriptTypes( [ScriptType.SCRIPTLET, ScriptType.GROOVY] )
class LocalExecutor extends Executor {

    /**
     * Environment variable to run the task wrapper script inside the Docker container
     * also on Linux and macOS hosts. On Windows this is always the case.
     */
    static final String DOCKER_WRAPPER_ENV = 'NXF_DOCKER_WRAPPER_IN_CONTAINER'

    @Override
    protected TaskMonitor createTaskMonitor() {
        return LocalPollingMonitor.create(session, config, name)
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
        return isFusionEnabled() || isDockerWrapperEnabled()
    }

    /**
     * When enabled, the whole task wrapper script is run inside the Docker container,
     * instead of running it on the host and launching the container from it.
     * This is always the case on Windows, since the host has no Bash.
     *
     * @return {@code true} when the task wrapper must run inside the Docker container
     */
    @Memoized
    boolean isDockerWrapperEnabled() {
        if( isFusionEnabled() )
            return false
        final config = session.getContainerConfig()
        if( !config?.isEnabled() || config.getEngine() != 'docker' )
            return false
        return isWindowsHost() || SysEnv.getBool(DOCKER_WRAPPER_ENV, false)
    }

    protected boolean isWindowsHost() {
        return SysHelper.isWindows()
    }

    @Override
    boolean isFusionEnabled() {
        return FusionHelper.isFusionEnabled(session)
    }
}

