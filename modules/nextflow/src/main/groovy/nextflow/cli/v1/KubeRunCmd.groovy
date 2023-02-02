/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.cli.v1

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.cli.KubeRunImpl

/**
 * Extends `run` command to support Kubernetes deployment
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = 'Execute a workflow in a Kubernetes cluster (experimental)')
class KubeRunCmd extends RunCmd implements KubeRunImpl.Options {

    @Parameter(names = ['-head-cpus'], description = 'Specify number of CPUs requested for the Nextflow driver pod')
    int headCpus

    @Parameter(names = ['-head-image'], description = 'Specify the container image for the Nextflow driver pod')
    String headImage

    @Parameter(names = ['-head-memory'], description = 'Specify amount of memory requested for the Nextflow driver pod')
    String headMemory

    @Parameter(names = ['-head-prescript'], description = 'Specify script to be run before nextflow run starts')
    String headPreScript

    @Parameter(names = ['-n','-namespace'], description = 'Specify the K8s namespace to use')
    String namespace

    @Parameter(names = ['-pod-image'], description = 'Alias for -head-image (deprecated)')
    String podImage

    @Parameter(names = ['-remoteConfig'], description = 'Add the specified file from the K8s cluster to configuration set', hidden = true)
    List<String> remoteConfig

    @Parameter(names = ['-remoteProfile'], description = 'Choose a configuration profile in the remoteConfig')
    String remoteProfile

    @Parameter(names = ['-v','-volume-mount'], description = 'Volume claim mounts eg. my-pvc:/mnt/path')
    List<String> volumeMounts

    @Override
    String getName() { 'kuberun' }

    @Override
    void run() {
        new KubeRunImpl(this).run()
    }

}
