/*
 * Copyright 2023, Seqera Labs
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

package nextflow.cli.v2

import groovy.transform.CompileStatic
import nextflow.cli.KubeRunImpl
import picocli.CommandLine.Command
import picocli.CommandLine.Option

/**
 * CLI `kuberun` sub-command (v2)
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Command(
    name = 'kuberun',
    description = 'Execute a workflow in a Kubernetes cluster (experimental)'
)
class KubeRunCmd extends RunCmd implements KubeRunImpl.Options {

    @Option(names = ['--head-cpus'], description = 'Specify number of CPUs requested for the Nextflow driver pod')
    int headCpus

    @Option(names = ['--head-image'], description = 'Specify the container image for the Nextflow driver pod')
    String headImage

    @Option(names = ['--head-memory'], description = 'Specify amount of memory requested for the Nextflow driver pod')
    String headMemory

    @Option(names = ['--head-prescript'], description = 'Specify script to be run before nextflow run starts')
    String headPreScript

    @Option(names = ['-n','--namespace'], description = 'Specify the K8s namespace to use')
    String namespace

    @Option(names = ['--pod-image'], description = 'Alias for -head-image (deprecated)')
    String podImage

    @Option(names = ['-v','--volume-mount'], description = 'Volume claim mounts eg. my-pvc:/mnt/path')
    List<String> volumeMounts

    @Option(names = [ '--remote-config'], description = 'Add the specified file from the K8s cluster to configuration set', hidden = true )
    List<String> remoteConfig

    @Option(names = ['--remote-profile'], description = 'Choose a configuration profile in the remoteConfig')
    String remoteProfile

    @Override
    Integer call() {
        new KubeRunImpl(this).run()
        return 0
    }

}
