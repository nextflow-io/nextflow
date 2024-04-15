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
 */

package nextflow.ga4gh.tes.executor

import groovy.transform.CompileStatic
import nextflow.executor.BashWrapperBuilder
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun

/**
 * Bash builder adapter to manage TES specific tasks
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TesBashBuilder extends BashWrapperBuilder {

    TesBashBuilder(TaskRun task, String remoteBinDir) {
        super(new TaskBean(task), new TesFileCopyStrategy(remoteBinDir))
    }

    TesBashBuilder(TaskBean task, String remoteBinDir) {
        super(task, new TesFileCopyStrategy(remoteBinDir))
    }

}
