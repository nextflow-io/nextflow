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
package nextflow.processor.hash

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
/**
 * Task hasher strategy that hashes the project bin scripts referenced in the
 * task script, on top of the V2 hashing behaviour.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskHasherV3 extends TaskHasherV2 {

    TaskHasherV3(TaskRun task) {
        super(task)
    }

    @Override
    protected List<Object> collectKeys() {
        final keys = super.collectKeys()
        // add bin scripts referenced in the task script
        final binEntries = getTaskBinEntries(task.source)
        if( binEntries ) {
            log.trace "Task: ${task.processor.name} > Adding scripts on project bin path: ${-> binEntries.join('; ')}"
            keys.addAll(binEntries)
        }
        return keys
    }

    @Memoized
    protected List<Path> getTaskBinEntries(String script) {
        return task.getTaskBinEntries(script)
    }
}
