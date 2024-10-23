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

package nextflow.fusion

import groovy.transform.CompileStatic
import nextflow.executor.Executor
import nextflow.processor.TaskRun
/**
 * Implements commons logic to handle fusion based tasks.
 * This trait is expected to be used by a sub-class of {@link nextflow.processor.TaskHandler}
 *
 * See {@link nextflow.executor.local.LocalTaskHandler}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
trait FusionAwareTask {

    abstract TaskRun getTask()

    private Executor getExecutor0() { getTask().getProcessor().getExecutor() }

    private Boolean fusionEnabled

    private FusionScriptLauncher fusionLauncher


    boolean fusionEnabled() {
        if( fusionEnabled==null ) {
            fusionEnabled = getExecutor0().isFusionEnabled()
        }
        return fusionEnabled
    }

    FusionConfig fusionConfig() {
        return FusionConfig.getConfig()
    }

    FusionScriptLauncher fusionLauncher() {
        if( fusionLauncher==null ) {
            fusionLauncher = fusionEnabled()
                    ? FusionScriptLauncher.create(task.toTaskBean(), task.workDir.scheme)
                    : null
        }
        return fusionLauncher
    }

    List<String> fusionSubmitCli() {
        return fusionLauncher().fusionSubmitCli(task)
    }

}
