/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.k8s.K8sDriverLauncher
/**
 * Extends `run` command to support Kubernetes deployment
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Execute a workflow in a Kubernetes cluster (experimental)")
class CmdKubeRun extends CmdRun {

    static private String POD_NAME = /[a-z0-9]([-a-z0-9]*[a-z0-9])?(\.[a-z0-9]([-a-z0-9]*[a-z0-9])?)*/

    /**
     * One or more volume claims mounts
     */
    @Parameter(names = ['-v','-volume-mount'], description = 'Volume claim mounts eg. my-pvc:/mnt/path')
    List<String> volMounts

    @Parameter(names = ['-n','-namespace'], description = 'Specify the K8s namespace to use')
    String namespace

    @Parameter(names = '-pod-image', description = 'Specify the container image for the Nextflow pod')
    String podImage

    @Override
    String getName() { 'kuberun' }

    @Override
    protected void checkRunName() {
        if( runName && !runName.matches(POD_NAME) )
            throw new AbortOperationException("Not a valid K8s pod name -- It can only contain lower case alphanumeric characters, '-' or '.', and must start and end with an alphanumeric character")
        super.checkRunName()
        runName = runName.replace('_','-')
    }

    @Override
    void run() {
        final scriptArgs = (args?.size()>1 ? args[1..-1] : []) as List<String>
        final pipeline = stdin ? '-' : ( args ? args[0] : null )
        if( !pipeline )
            throw new AbortOperationException("No project name was specified")

        checkRunName()
        new K8sDriverLauncher(cmd: this, runName: runName, podImage: podImage).run(pipeline, scriptArgs)
    }

}
