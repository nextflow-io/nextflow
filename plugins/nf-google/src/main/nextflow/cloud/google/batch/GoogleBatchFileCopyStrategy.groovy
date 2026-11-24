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

package nextflow.cloud.google.batch

import groovy.transform.CompileStatic
import nextflow.executor.SimpleFileCopyStrategy
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun

/**
 * File copy strategy for Google Batch tasks.
 * <p>
 * The launcher stages the task's own command files ({@code .command.sh},
 * {@code .command.run}, {@code .command.in}) as inputs, so that they exist in
 * the scratch directory. These are always copied, whatever the {@code stageInMode}:
 * a symlink back to the work dir would make a {@code path('.command.*')} output
 * be copied onto itself at unstage time (see #6213).
 * Every other input keeps the configured stage-in mode.
 */
@CompileStatic
class GoogleBatchFileCopyStrategy extends SimpleFileCopyStrategy {

    private static final Set<String> TASK_FILES = Set.of(TaskRun.CMD_SCRIPT, TaskRun.CMD_RUN, TaskRun.CMD_INFILE)

    GoogleBatchFileCopyStrategy(TaskBean bean) {
        super(bean)
    }

    @Override
    protected String stageInCommand(String source, String target, String mode) {
        return super.stageInCommand(source, target, TASK_FILES.contains(target) ? 'copy' : mode)
    }
}
