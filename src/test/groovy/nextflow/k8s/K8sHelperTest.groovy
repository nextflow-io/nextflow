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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sHelperTest extends Specification {

    def setup() {
        K8sHelper.VOLUMES.set(0)
    }

    def 'should get storage spec for volume claims' () {

        given:
        final REQ = [podName: 'foo', imageName: 'busybox', workDir: '/path', command: ['echo']]

        when:
        REQ.volumeClaims = [ first: [mountPath: '/work'], second: [mountPath: '/data']]
        def spec = K8sHelper.createPodSpec(REQ)
        then:
        spec ==  [ apiVersion: 'v1',
                   kind: 'Pod',
                   metadata: [name:'foo', namespace:'default'],
                   spec: [
                           restartPolicy:'Never',
                           containers:[
                                   [name:'foo',
                                    image:'busybox',
                                    command: ['echo'],
                                    workingDir:'/path',
                                    volumeMounts:[
                                            [name:'vol-1', mountPath:'/work'],
                                            [name:'vol-2', mountPath:'/data']] ]
                           ],
                           volumes:[
                                   [name:'vol-1', persistentVolumeClaim:[claimName:'first']],
                                   [name:'vol-2', persistentVolumeClaim:[claimName:'second']] ]
                   ]

        ]

    }

    def 'should get config map mounts' () {

        given:
        final REQ = [podName: 'foo', imageName: 'busybox', workDir: '/path', command: ['echo']]

        when:
        REQ.configMounts = [cfg1: '/etc/config', data2: '/data/path']

        def spec = K8sHelper.createPodSpec(REQ)
        then:
        spec ==  [ apiVersion: 'v1',
                   kind: 'Pod',
                   metadata: [name:'foo', namespace:'default'],
                   spec: [
                           restartPolicy:'Never',
                           containers:[
                                   [name:'foo',
                                    image:'busybox',
                                    command: ['echo'],
                                    workingDir:'/path',
                                    volumeMounts:[
                                            [name:'vol-1', mountPath:'/etc/config'],
                                            [name:'vol-2', mountPath:'/data/path']] ]
                           ],
                           volumes:[
                                   [name:'vol-1', configMap:[name:'cfg1']],
                                   [name:'vol-2', configMap:[name:'data2']] ]
                   ]

        ]

    }


    def 'should get host path mounts' () {
        given:
        final REQ = [podName: 'foo', imageName: 'busybox', workDir: '/path', command: ['echo']]

        when:
        REQ.hostMounts = ['/tmp':'/scratch', '/host/data':'/mnt/container']
        def spec = K8sHelper.createPodSpec(REQ)
        then:
        spec ==  [ apiVersion: 'v1',
                   kind: 'Pod',
                   metadata: [name:'foo', namespace:'default'],
                   spec: [
                           restartPolicy:'Never',
                           containers:[
                                   [name:'foo',
                                    image:'busybox',
                                    command: ['echo'],
                                    workingDir:'/path',
                                    volumeMounts:[
                                            [name:'vol-1', mountPath:'/scratch'],
                                            [name:'vol-2', mountPath:'/mnt/container']] ]
                           ],
                           volumes:[
                                   [name:'vol-1', hostPath: [path:'/tmp']],
                                   [name:'vol-2', hostPath: [path:'/host/data']] ]
                   ]

        ]

    }

    def 'should create pod spec' () {

        when:
        def spec = K8sHelper.createPodSpec([
                        podName: 'foo',
                        imageName: 'busybox',
                        command: ['echo', 'hello'],
                        workDir: '/some/work/dir'])
        then:
        spec ==  [ apiVersion: 'v1',
                   kind: 'Pod',
                   metadata: [name:'foo', namespace:'default'],
                   spec: [
                           restartPolicy:'Never',
                           containers:[
                                   [name:'foo',
                                    image:'busybox',
                                    command:['echo', 'hello'],
                                    workingDir:'/some/work/dir'
                                   ]
                           ]
                   ]
                ]

    }


    def 'should set namespace and labels' () {

        when:
        def spec = K8sHelper.createPodSpec([
                podName: 'foo',
                imageName: 'busybox',
                command: ['sh', '-c', 'echo hello'],
                workDir: '/some/work/dir',
                namespace: 'xyz',
                labels: [app: 'myApp', runName: 'something']
        ])
        then:
        spec ==  [ apiVersion: 'v1',
                   kind: 'Pod',
                   metadata: [name:'foo', namespace:'xyz', labels:[app: 'myApp', runName: 'something']],
                   spec: [
                           restartPolicy:'Never',
                           containers:[
                                   [name:'foo',
                                    image:'busybox',
                                    command: ['sh', '-c', 'echo hello'],
                                    workingDir:'/some/work/dir'
                                   ]
                           ]
                   ]
        ]
    }

    def 'should set resources and env' () {

        when:
        def spec = K8sHelper.createPodSpec([
                podName: 'foo',
                imageName: 'busybox',
                command: 'echo hello',
                workDir: '/some/work/dir',
                env: [ALPHA: 'hello', DELTA: 'world'],
                cpus: 8,
                memory: '100Gi'
        ])
        then:
        spec ==  [ apiVersion: 'v1',
                   kind: 'Pod',
                   metadata: [name:'foo', namespace:'default'],
                   spec: [
                           restartPolicy:'Never',
                           containers:[
                                   [name:'foo',
                                    image:'busybox',
                                    command:['/bin/bash', '-c', 'echo hello'],
                                    workingDir:'/some/work/dir',
                                    env: [
                                      [name:'ALPHA', value:'hello'],
                                      [name:'DELTA', value:'world']
                                    ],
                                    resources:[limits:[cpu:8, memory:'100Gi'] ]
                                   ]
                           ]
                   ]
        ]
    }
}
