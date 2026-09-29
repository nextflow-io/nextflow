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
 * paths the filesystem tool may read — the work dir, the SOURCES of the task's staged
 * {@code Path} inputs (the guard resolves the stage-in symlink, so the source is what the
 * containment test actually sees), and the outputs of modules run during this invocation.
 * Created per input record by the agent operator and
 * threaded to {@link ModuleToolBridge} via a ThreadLocal, so the shared,
 * pre-ignition bridge holds no per-record state.
 *
 * <p>An entry may be a FILE or a directory: containment is by path prefix, so a file
 * entry grants exactly that file while a directory entry grants its whole subtree.
 * Module outputs are added as files precisely so that reading one does not also grant
 * its siblings — see {@link ModuleToolBridge#collectPathsFromValue}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DispatchContext {
    final Path workDir
    final Set<Path> readablePaths

    DispatchContext(Path workDir) {
        this.workDir = workDir
        this.readablePaths = ConcurrentHashMap.newKeySet()
        if( workDir != null )
            this.readablePaths.add(workDir)
    }

    void addReadablePath(Path path) {
        if( path != null )
            readablePaths.add(path)
    }
}
