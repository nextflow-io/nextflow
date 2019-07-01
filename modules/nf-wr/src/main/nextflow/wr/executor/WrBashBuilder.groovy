/*
 * Copyright 2019, Genome Research Limited
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

package nextflow.wr.executor

import groovy.transform.CompileStatic
import java.nio.file.Path

import nextflow.executor.BashWrapperBuilder
import nextflow.processor.TaskRun
import nextflow.processor.TaskBean

/**
 * Bash builder adapter to manage wr specific tasks
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * based on TesBashBuilder by Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class WrBashBuilder extends BashWrapperBuilder {

    WrBashBuilder(TaskRun task, Path remoteBinDir) {
        super(new TaskBean(task), new WrFileCopyStrategy(new TaskBean(task), remoteBinDir))
    }

    WrBashBuilder(TaskBean bean, Path remoteBinDir) {
        super(bean, new WrFileCopyStrategy(bean, remoteBinDir))
    }

    protected boolean alwaysTryToUnstage() {
        return true
    }

}