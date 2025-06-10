/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.k8s.cli

import groovy.transform.CompileStatic
import nextflow.cli.CmdKubeRun
import nextflow.k8s.K8sDriverLauncher

/**
 * Kuberun command implementation logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class KubeCommandImpl implements CmdKubeRun.KubeCommand {

    @Override
    int run(CmdKubeRun cmd, String pipeline, List<String> args) {
        // create
        final driver = new K8sDriverLauncher(
            cmd: cmd,
            runName: cmd.runName,
            headImage: cmd.headImage,
            background: cmd.background(),
            headCpus: cmd.headCpus,
            headMemory: cmd.headMemory,
            headPreScript: cmd.headPreScript,
            plugins: cmd.plugins )
        // run it
        driver.run(pipeline, args)
        // return exit code
        return driver.shutdown()
    }
    
}
