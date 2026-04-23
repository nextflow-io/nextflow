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
package nextflow.executor

import groovy.transform.CompileStatic
import nextflow.executor.local.LocalExecutor
import nextflow.script.ScriptType
import nextflow.util.ServiceName

/**
 * Executes tasks locally inside Apple containers (Apple Silicon Macs).
 *
 * Equivalent to the {@code local} executor but automatically selects the
 * {@code apple} container engine, so users only need to set
 * {@code process.executor = 'apple_container'} and supply a container image.
 *
 * @author Joon-Klaps
 */
@CompileStatic
@ServiceName('apple_container')
@SupportedScriptTypes([ScriptType.SCRIPTLET])
class AppleContainerExecutor extends LocalExecutor {

    @Override
    String containerConfigEngine() {
        return 'apple'
    }

}
