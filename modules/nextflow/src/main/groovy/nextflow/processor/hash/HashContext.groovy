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
import nextflow.Session
import nextflow.processor.TaskHasher
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun

/**
 * Everything a Contributor may read while extracting values for a task.
 *
 * The two derived helpers are delegated to the existing TaskHasher rather than
 * duplicated, so the interpreter and the byte-exactness oracle cannot diverge on
 * them. The helper is injectable so a test can stub both sides identically.
 */
@CompileStatic
class HashContext {

    final TaskRun task

    final TaskProcessor processor

    final Session session

    /** Services injected by a plugin-supplied spec, keyed by name. */
    final Map<String,Object> services

    private final TaskHasher helper

    HashContext(TaskRun task) {
        this(task, new TaskHasher(task), [:])
    }

    HashContext(TaskRun task, TaskHasher helper) {
        this(task, helper, [:])
    }

    HashContext(TaskRun task, TaskHasher helper, Map<String,Object> services) {
        this.task = task
        this.processor = task.processor
        this.session = task.processor.session
        this.helper = helper
        this.services = services
    }

    Map<String,Object> globalVars() {
        return helper.getTaskGlobalVars()
    }

    List<Path> binEntries() {
        return helper.getTaskBinEntries(task.source)
    }
}
