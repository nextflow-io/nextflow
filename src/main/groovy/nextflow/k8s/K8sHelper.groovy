/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.k8s

import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.PackageScope

/**
 * Helper class to handle Kubernetes requests
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sHelper {

    static @PackageScope AtomicInteger VOLUMES = new AtomicInteger()

    /**
     * Create a pod spec
     *
     * @param params Accepted parameters:
     *      podName
     *      imageName
     *      command
     *      workDir
     *      namespace
     *      cpus (int)
     *      env (Map of key-value pairs)
     *      memory (string using k8s format)
     *      labels (Map of key-value pairs)
     *      volumeClaims (Map of volume claims)
     *      configMounts (Map of config mounts)
     *      hostMounts (Map of hostPath mounts)
     *      serviceAccount
     * @return
     */
    static Map createPodSpec( Map params ) {
        assert params.podName, 'Missing K8s podName parameter'
        assert params.imageName, 'Missing K8s imageName parameter'
        assert params.command instanceof List || params.command instanceof CharSequence, "Missing or invalid K8s command parameter: $params.command"

        final namespace = params.namespace ?: 'default'
        final restart = params.restart ?: 'Never'

        final labels = params.labels instanceof Map ? params.labels : [:]
        final cmd = params.command instanceof List ? params.command : ['/bin/bash','-c', params.command.toString()]
        final env = []
        if( params.env instanceof Map ) {
            params.env.each { k, v -> env << [name: k, value: v] }
        }

        final res = [:]
        if( params.cpus )
            res.cpu = params.cpus
        if( params.memory )
            res.memory = params.memory

        final container = [
                name: params.podName,
                image: params.imageName,
                command: cmd
        ]
        if( params.workDir )
            container.workingDir = params.workDir.toString()

        final spec = [
                restartPolicy: restart,
                containers: [ container ],
        ]

        if( params.serviceAccount )
            spec.serviceAccountName = params.serviceAccount

        final pod = [
                apiVersion: 'v1',
                kind: 'Pod',
                metadata: [
                        name: params.podName,
                        namespace: namespace
                ],
                spec: spec
        ]

        // add environment
        if( env )
            container.env = env

        // add resources
        if( res ) {
            container.resources = [limits: res]
        }

        // add labels
        if( labels )
            pod.metadata.labels = labels

        // add storage definitions ie. volumes and mounts
        final mounts = []
        final volumes = []

        // -- volume claims
        if( params.volumeClaims instanceof Map ) {
            params.volumeClaims.each { String claimName, Map entry ->
                final name = nextVolName()
                final claim = [name: name, mountPath: entry.mountPath ]
                mounts << claim
                volumes << [name: name, persistentVolumeClaim: [claimName: claimName]]
            }
        }

        // -- configmap volumes
        if( params.configMounts instanceof Map ) {
            params.configMounts.each { cfg, path ->
                final name = nextVolName()
                mounts << [name: name, mountPath: path]
                volumes << [name: name, configMap: [name: cfg]]
            }
        }

        if( params.hostMounts instanceof Map ) {
            params.hostMounts.each { String hostPath, String mountPath ->
                final name = nextVolName()
                mounts << [name: name, mountPath: mountPath]
                volumes << [name: name, hostPath: [path: hostPath]]
            }
        }

        if( volumes )
            pod.spec.volumes = volumes
        if( mounts )
            container.volumeMounts = mounts

        return pod
    }


    /**
     * @return A sequential volume unique identifier
     */
    static protected String nextVolName() {
        "vol-${VOLUMES.incrementAndGet()}".toString()
    }

}
