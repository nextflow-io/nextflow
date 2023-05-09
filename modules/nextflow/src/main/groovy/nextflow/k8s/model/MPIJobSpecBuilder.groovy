/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.k8s.model

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.util.MemoryUnit

/**
 * Builder for a K8s MPIJob specification
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Slf4j
class MPIJobSpecBuilder extends PodSpecBuilder {

    Integer workers = 1

    String sshAuthMountPath

    boolean waitForWorkers = false

    @Override
    ResourceSpecBuilder withPodOptions(PodOptions opts) {
        super.withPodOptions(opts)
        withMpiOptions(opts.unmatched.mpi as Map)
    }

    ResourceSpecBuilder withMpiOptions(Map opts) {
        if( opts.workers )
            workers = opts.workers as Integer
        if( opts.sshAuthMountPath )
            sshAuthMountPath = opts.sshAuthMountPath as String
        if( opts.waitForWorkers )
            waitForWorkers = opts.waitForWorkers as Boolean
        return this
    }

    @Override
    Map build() {
        // build metadata
        final metadata = super.buildMetadata()

        // save command and args
        final command = this.command
        final args = this.args

        // build worker spec
        withCommand([])
        withArgs([])

        final workerSpec = super.buildSpec()

        // build launcher spec
        withCommand(command)
        withArgs(args)
        withCpus(1)
        withMemory(MemoryUnit.of('2 GB'))

        final launcherSpec = super.buildSpec()

        // build MPIJob spec
        final spec = [
            slotsPerWorker: 1,
            runPolicy: [
                cleanPodPolicy: 'All'
            ],
            sshAuthMountPath: this.sshAuthMountPath,
            waitForWorkers: this.waitForWorkers,
            mpiReplicaSpecs: [
                Launcher: [
                    replicas: 1,
                    template: [
                        spec: launcherSpec
                    ]
                ],
                Worker: [
                    replicas: this.workers,
                    template: [
                        spec: workerSpec
                    ]
                ]
            ]
        ]

        return [
            apiVersion: 'kubeflow.org/v2beta1',
            kind: 'MPIJob',
            metadata: metadata,
            spec: spec
        ]
    }

}
