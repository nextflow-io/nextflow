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

package nextflow.cli

import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic

/**
 * Stub for the removed `kuberun` command.
 */
@CompileStatic
@Parameters
class CmdKubeRun extends CmdRemoved {

    static final public String NAME = 'kuberun'

    @Override
    String getName() { NAME }

    @Override
    String getMessage() {
        "The `kuberun` command has been removed -- use `nextflow plugin nf-kuberun:run` from the nf-kuberun plugin, which uses double-dash options (e.g. `--namespace`, `--volume-mount`), or launch with Fusion instead"
    }
}
