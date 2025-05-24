/*
 * Copyright 2021, Microsoft Corp
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

package nextflow.cloud.azure.batch

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.Executor
import nextflow.fusion.FusionAwareTask
import nextflow.fusion.FusionConfig
import nextflow.fusion.FusionScriptLauncher
import nextflow.processor.TaskRun

/**
 * Adapter class that wraps a TaskRun to implement the FusionAwareTask interface
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AzFusionTaskWrapper implements FusionAwareTask {

    private TaskRun task
    private FusionScriptLauncher fusionLauncher
    private Boolean fusionEnabled

    AzFusionTaskWrapper(TaskRun task) {
        this.task = task
    }

    @Override
    TaskRun getTask() {
        return task
    }

    @Override
    boolean fusionEnabled() {
        return true  // Always true since we only create this wrapper when fusion is enabled
    }

    @Override
    FusionConfig fusionConfig() {
        return FusionConfig.getConfig()
    }

    @Override
    FusionScriptLauncher fusionLauncher() {
        if (fusionLauncher == null) {
            fusionLauncher = FusionScriptLauncher.create(task.toTaskBean(), task.workDir.scheme)
        }
        return fusionLauncher
    }

    @Override
    List<String> fusionSubmitCli() {
        return fusionLauncher().fusionSubmitCli(task)
    }
} 