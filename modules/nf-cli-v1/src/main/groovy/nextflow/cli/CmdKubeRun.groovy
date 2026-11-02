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

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.exception.AbortOperationException

/**
 * Stub for the removed `kuberun` command, kept to report a clear error.
 * Omitting the command description hides it from the command listing.
 */
@CompileStatic
@Parameters
class CmdKubeRun extends CmdBase implements UsageAware {

    static final public String NAME = 'kuberun'

    @Parameter
    List<String> args

    String getName() { NAME }

    void run() {
        abort()
    }

    @Override
    void usage() {
        abort()
    }

    @Override
    void usage(List<String> args) {
        abort()
    }

    private static void abort() {
        throw new AbortOperationException("The `kuberun` command has been removed -- use the `nf-kuberun` plugin or launch with Fusion instead")
    }
}
