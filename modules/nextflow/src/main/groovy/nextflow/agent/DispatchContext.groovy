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
package nextflow.agent

import java.nio.file.Path
import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic

/**
 * Per-agent-invocation dispatch context: the sandbox work dir and the set of
 * directories the filesystem tool may read (the work dir plus the output dirs of
 * modules run during this invocation). Created per input record by the agent
 * operator and threaded to {@link ModuleToolBridge} via a ThreadLocal, so the
 * shared, pre-ignition bridge holds no per-record state.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DispatchContext {
    final Path workDir
    final Set<Path> readableDirs

    DispatchContext(Path workDir) {
        this.workDir = workDir
        this.readableDirs = ConcurrentHashMap.newKeySet()
        if( workDir != null )
            this.readableDirs.add(workDir)
    }

    void addReadableDir(Path dir) {
        if( dir != null )
            readableDirs.add(dir)
    }
}
