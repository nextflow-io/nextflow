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

package nextflow.k8s

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.Session
import nextflow.exception.NodeTerminationException
import nextflow.executor.fusion.FusionScriptLauncher
import nextflow.file.http.XPath
import nextflow.k8s.client.ClientConfig
import nextflow.k8s.client.K8sClient
import nextflow.k8s.client.K8sResponseException
import nextflow.k8s.client.K8sResponseJson
import nextflow.k8s.model.PodEnv
import nextflow.k8s.model.PodMountConfig
import nextflow.k8s.model.PodMountSecret
import nextflow.k8s.model.PodOptions
import nextflow.k8s.model.PodSpecBuilder
import nextflow.k8s.model.PodVolumeClaim
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.util.MemoryUnit
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sTaskHandlerTest extends Specification {

    def setup() {
        PodSpecBuilder.VOLUMES.set(0)
    }

    def 'should return a new pod with args' () {
        given:
        def WORK_DIR = Paths.get('/some/work/dir')
        def config = Mock(TaskConfig)
        def task = Mock(TaskRun)
        def client = Mock(K8sClient)
        def builder = Mock(K8sWrapperBuilder)
        def handler = Spy(new K8sTaskHandler(builder:builder, client: client))
        Map result

        when:
        result = handler.newSubmitRequest(task)
        then:
        _ * handler.fusionEnabled() >> false
        1 * handler.fixOwnership() >> false
        1 * handler.entrypointOverride() >> false
        1 * handler.getPodOptions() >> new PodOptions()
        1 * handler.getSyntheticPodName(task) >> 'nf-123'
        1 * handler.getLabels(task) >> [:]
        1 * handler.getAnnotations() >> [:]
        1 * handler.getContainerMounts() >> []
        1 * task.getContainer() >> 'debian:latest'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> config
        1 * config.getCpus() >> 0
        1 * config.getMemory() >> null
        1 * client.getConfig() >> new ClientConfig()
        result == [ apiVersion: 'v1',
                    kind: 'Pod',
                    metadata: [
                            name:'nf-123',
                            namespace:'default'
                    ],
                    spec: [
                            restartPolicy:'Never',
                            containers:[
                                    [name:'nf-123',
                                     image:'debian:latest',
                                     args:['/bin/bash', '-ue','/some/work/dir/.command.run'] ]
                            ]
                    ]
        ]

        when:
        result = handler.newSubmitRequest(task)
        then:
        _ * handler.fusionEnabled() >> false
        1 * handler.entrypointOverride() >> true
        1 * handler.getSyntheticPodName(task) >> 'nf-foo'
        1 * handler.getLabels(task) >> [sessionId:'xxx']
        1 * handler.getAnnotations() >>  [evict: 'false']
        1 * handler.getPodOptions() >> new PodOptions()
        1 * handler.getContainerMounts() >> []
        1 * handler.fixOwnership() >> true
        1 * handler.getOwner() >> '501:502'
        1 * task.getContainer() >> 'debian:latest'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> config
        1 * config.getCpus() >> 1
        1 * config.getMemory() >> null
        1 * client.getConfig() >> new ClientConfig()
        result == [ apiVersion: 'v1',
                    kind: 'Pod',
                    metadata: [name:'nf-foo', namespace:'default', labels: [sessionId: 'xxx'], annotations: [evict: 'false']],
                    spec: [
                            restartPolicy:'Never',
                            containers:[
                                    [name:'nf-foo',
                                     image:'debian:latest',
                                     command:['/bin/bash', '-ue','/some/work/dir/.command.run'],
                                     resources:[ requests: [cpu:1] ],
                                     env: [  [name:'NXF_OWNER', value:'501:502'] ]
                                    ]
                            ]
                    ]
        ]


        when:
        result = handler.newSubmitRequest(task)
        then:
        _ * handler.fusionEnabled() >> false
        1 * handler.fixOwnership() >> false
        1 * handler.entrypointOverride() >> true
        1 * handler.getSyntheticPodName(task) >> 'nf-abc'
        1 * handler.getLabels(task) >> [:]
        1 * handler.getAnnotations() >> [:]
        1 * handler.getPodOptions() >> new PodOptions()
        1 * handler.getContainerMounts() >> []
        1 * task.getContainer() >> 'user/alpine:1.0'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> config
        1 * config.getCpus() >> 4
        1 * config.getMemory() >> MemoryUnit.of('16GB')
        1 * client.getConfig() >> new ClientConfig(namespace: 'namespace-x')
        result == [ apiVersion: 'v1',
                    kind: 'Pod',
                    metadata: [name:'nf-abc', namespace:'namespace-x' ],
                    spec: [
                            restartPolicy:'Never',
                            containers:[
                                    [name:'nf-abc',
                                     image:'user/alpine:1.0',
                                     command:['/bin/bash', '-ue', '/some/work/dir/.command.run'],
                                     resources:[ requests: [cpu:4, memory:'16384Mi'], limits:[memory:'16384Mi'] ]
                                    ]
                            ]
                    ]
        ]

    }

    def 'should return a new pod request with no storage' () {
        given:
        def WORK_DIR = Paths.get('/some/work/dir')
        def config = Mock(TaskConfig)
        def task = Mock(TaskRun)
        def client = Mock(K8sClient)
        def builder = Mock(K8sWrapperBuilder)
        def handler = Spy(new K8sTaskHandler(builder: builder, client:client))
        Map result

        when:
        result = handler.newSubmitRequest(task)
        then:
        _ * handler.fusionEnabled() >> false
        1 * handler.fixOwnership() >> false
        1 * handler.entrypointOverride() >> true
        1 * handler.getPodOptions() >> new PodOptions()
        1 * handler.getSyntheticPodName(task) >> 'nf-123'
        1 * handler.getLabels(task) >> [foo: 'bar', hello: 'world']
        1 * handler.getAnnotations() >> [fooz: 'barz', ciao: 'mondo']
        1 * handler.getContainerMounts() >> []
        1 * task.getContainer() >> 'debian:latest'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> config
        1 * config.getCpus() >> 0
        1 * config.getMemory() >> null
        1 * client.getConfig() >> new ClientConfig()
        result == [ apiVersion: 'v1',
                    kind: 'Pod',
                    metadata: [
                            name:'nf-123',
                            namespace:'default',
                            labels:[ foo:'bar', hello: 'world'],
                            annotations:[ fooz:'barz', ciao: 'mondo']
                    ],
                    spec: [
                            restartPolicy:'Never',
                            containers:[
                                    [name:'nf-123',
                                     image:'debian:latest',
                                     command:['/bin/bash', '-ue','/some/work/dir/.command.run'] ]
                            ]
                    ]
                ]

        when:
        result = handler.newSubmitRequest(task)
        then:
        _ * handler.fusionEnabled() >> false
        1 * handler.entrypointOverride() >> true
        1 * handler.getSyntheticPodName(task) >> 'nf-foo'
        1 * handler.getLabels(task) >> [sessionId:'xxx']
        1 * handler.getAnnotations() >>  [evict: 'false']
        1 * handler.getPodOptions() >> new PodOptions()
        1 * handler.getContainerMounts() >> []
        1 * handler.fixOwnership() >> true
        1 * handler.getOwner() >> '501:502'
        1 * task.getContainer() >> 'debian:latest'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> config
        1 * config.getCpus() >> 1
        1 * config.getMemory() >> null
        1 * client.getConfig() >> new ClientConfig()
        result == [ apiVersion: 'v1',
                    kind: 'Pod',
                    metadata: [name:'nf-foo', namespace:'default', labels: [sessionId: 'xxx'], annotations: [evict: 'false']],
                    spec: [
                            restartPolicy:'Never',
                            containers:[
                                    [name:'nf-foo',
                                     image:'debian:latest',
                                     command:['/bin/bash', '-ue','/some/work/dir/.command.run'],
                                     resources:[ requests: [cpu:1] ],
                                     env: [  [name:'NXF_OWNER', value:'501:502'] ]
                                    ]
                            ]
                    ]
        ]

        when:
        result = handler.newSubmitRequest(task)
        then:
        _ * handler.fusionEnabled() >> false
        1 * handler.fixOwnership() >> false
        1 * handler.entrypointOverride() >> true
        1 * handler.getSyntheticPodName(task) >> 'nf-abc'
        1 * handler.getLabels(task) >> [:]
        1 * handler.getAnnotations() >> [:]
        1 * handler.getPodOptions() >> new PodOptions()
        1 * handler.getContainerMounts() >> []
        1 * task.getContainer() >> 'user/alpine:1.0'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> config
        1 * config.getCpus() >> 4
        1 * config.getMemory() >> MemoryUnit.of('16GB')
        1 * client.getConfig() >> new ClientConfig(namespace: 'namespace-x')
        result == [ apiVersion: 'v1',
                    kind: 'Pod',
                    metadata: [name:'nf-abc', namespace:'namespace-x' ],
                    spec: [
                            restartPolicy:'Never',
                            containers:[
                                    [name:'nf-abc',
                                     image:'user/alpine:1.0',
                                     command:['/bin/bash', '-ue', '/some/work/dir/.command.run'],
                                     resources:[ requests: [cpu:4, memory:'16384Mi'], limits: [memory:'16384Mi'] ]
                                    ]
                            ]
                    ]
        ]

    }

    def 'should create a pod with specified client configs' () {

        given:
        def WORK_DIR = Paths.get('/some/work/dir')
        def task = Mock(TaskRun)
        def client = Mock(K8sClient)
        def builder = Mock(K8sWrapperBuilder)
        def config = Mock(ClientConfig)
        def handler = Spy(new K8sTaskHandler(builder: builder, client: client))
        Map result

        when:
        result = handler.newSubmitRequest(task)
        then:
        _ * handler.fusionEnabled() >> false
        1 * handler.fixOwnership() >> false
        1 * handler.entrypointOverride() >> true
        1 * handler.getSyntheticPodName(task) >> 'nf-123'
        1 * handler.getPodOptions() >> new PodOptions()
        1 * handler.getLabels(task) >> [:]
        1 * handler.getAnnotations() >> [:]
        1 * handler.getContainerMounts() >> []
        1 * task.getContainer() >> 'debian:latest'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> new TaskConfig()
        1 * client.getConfig() >> config
        1 * config.getNamespace() >> 'just-a-namespace'
        1 * config.getServiceAccount() >> 'pedantic-kallisto'

        result == [ apiVersion: 'v1',
                    kind: 'Pod',
                    metadata: [name:'nf-123', namespace:'just-a-namespace' ],
                    spec: [
                            serviceAccountName: 'pedantic-kallisto',
                            restartPolicy:'Never',
                            containers:[
                                    [name:'nf-123',
                                     image:'debian:latest',
                                     command:['/bin/bash', '-ue','/some/work/dir/.command.run'],
                                     resources:[requests:[cpu:1]]
                                    ]
                            ]
                    ]
        ]

    }

    def 'should create a pod with specified pod options' () {

        given:
        def WORK_DIR = Paths.get('/some/work/dir')
        def task = Mock(TaskRun)
        def client = Mock(K8sClient)
        def builder = Mock(K8sWrapperBuilder)
        def config = Mock(TaskConfig)
        def handler = Spy(new K8sTaskHandler(builder:builder, client:client))
        def podOptions = Mock(PodOptions)
        and:
        podOptions.automountServiceAccountToken >> true
        Map result

        when:
        result = handler.newSubmitRequest(task)
        then:
        1 * client.getConfig() >> new ClientConfig()
        _ * handler.fusionEnabled() >> false
        1 * handler.fixOwnership() >> false
        1 * handler.entrypointOverride() >> true
        1 * handler.getSyntheticPodName(task) >> 'nf-123'
        1 * handler.getLabels(task) >> [:]
        1 * handler.getAnnotations() >> [:]
        1 * handler.getPodOptions() >> podOptions
        1 * handler.getContainerMounts() >> []
        1 * task.getContainer() >> 'debian:latest'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> config
        2 * podOptions.getEnvVars() >> [ PodEnv.value('FOO','bar') ]
        2 * podOptions.getMountSecrets() >> [ new PodMountSecret('my-secret/key-z', '/data/secret.txt') ]
        2 * podOptions.getMountConfigMaps() >> [ new PodMountConfig('my-data/key-x', '/etc/file.txt') ]

        result == [
            apiVersion: 'v1',
            kind: 'Pod',
            metadata: [name:'nf-123', namespace:'default' ],
            spec: [
                restartPolicy: 'Never',
                containers: [
                    [
                        name: 'nf-123',
                        image: 'debian:latest',
                        command: ['/bin/bash', '-ue','/some/work/dir/.command.run'],
                        env: [[name:'FOO', value:'bar']],
                        volumeMounts: [
                            [name:'vol-1', mountPath:'/etc'],
                            [name:'vol-2', mountPath:'/data']
                        ]
                    ]
                ],
                volumes:[
                    [name:'vol-1', configMap:[name:'my-data', items:[[key:'key-x', path:'file.txt']]]],
                    [name:'vol-2', secret:[secretName:'my-secret', items:[[key:'key-z', path:'secret.txt']]]]
                ]
            ]
        ]

    }

    def 'should create a request with vols and mounts' () {

        given:
        def WORK_DIR = Paths.get('/some/work/dir')
        def config = Mock(TaskConfig)
        def task = Mock(TaskRun)
        def client = Mock(K8sClient)
        def builder = Mock(K8sWrapperBuilder)
        def handler = Spy(new K8sTaskHandler(builder:builder, client:client))
        and:
        Map result

        def podOptions = Mock(PodOptions)
        def CLAIMS = [ new PodVolumeClaim('first','/work'), new PodVolumeClaim('second','/data') ]
        podOptions.automountServiceAccountToken >> true

        when:
        result = handler.newSubmitRequest(task)
        then:
        _ * handler.fusionEnabled() >> false
        1 * handler.fixOwnership() >> false
        1 * handler.entrypointOverride() >> true
        1 * handler.getSyntheticPodName(task) >> 'nf-123'
        1 * handler.getContainerMounts() >> []
        1 * handler.getLabels(task) >> [:]
        1 * handler.getAnnotations() >> [:]
        1 * handler.getPodOptions() >> podOptions
        1 * task.getContainer() >> 'debian:latest'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> config
        1 * config.getCpus() >> 0
        1 * config.getMemory() >> null
        1 * client.getConfig() >> new ClientConfig()
        2 * podOptions.getVolumeClaims() >> CLAIMS

        result == [
            apiVersion: 'v1',
            kind: 'Pod',
            metadata: [name:'nf-123', namespace:'default'],
            spec: [
                restartPolicy: 'Never',
                containers: [
                    [
                        name: 'nf-123',
                        image: 'debian:latest',
                        command: ['/bin/bash', '-ue', '/some/work/dir/.command.run'],
                        volumeMounts: [
                            [name:'vol-1', mountPath:'/work'],
                            [name:'vol-2', mountPath:'/data']
                        ]
                    ]
                ],
                volumes: [
                    [name:'vol-1', persistentVolumeClaim:[claimName: 'first']],
                    [name:'vol-2', persistentVolumeClaim:[claimName: 'second']]
                ]
            ]
        ]


        when:
        result = handler.newSubmitRequest(task)
        then:
        _ * handler.fusionEnabled() >> false
        1 * handler.fixOwnership() >> false
        1 * handler.entrypointOverride() >> true
        1 * handler.getSyntheticPodName(task) >> 'nf-123'
        1 * handler.getContainerMounts() >> ['/tmp', '/data']
        1 * handler.getLabels(task) >> [:]
        1 * handler.getAnnotations() >> [:]
        1 * handler.getPodOptions() >> new PodOptions()
        1 * task.getContainer() >> 'debian:latest'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> config
        1 * config.getCpus() >> 0
        1 * config.getMemory() >> null
        1 * client.getConfig() >> new ClientConfig()

        result == [
            apiVersion: 'v1',
            kind: 'Pod',
            metadata: [name:'nf-123', namespace:'default'],
            spec: [
                restartPolicy: 'Never',
                containers: [
                    [
                        name: 'nf-123',
                        image: 'debian:latest',
                        command: ['/bin/bash', '-ue', '/some/work/dir/.command.run'],
                        volumeMounts: [
                            [name:'vol-3', mountPath:'/tmp'],
                            [name:'vol-4', mountPath: '/data']
                        ]
                    ]
                ],
                volumes: [
                    [name:'vol-3', hostPath:[path:'/tmp']],
                    [name:'vol-4', hostPath:[path:'/data']]
                ]
            ]
        ]

    }

    def 'should submit a pod' () {

        given:
        def task = Mock(TaskRun)
        def client = Mock(K8sClient)
        def builder = Mock(K8sWrapperBuilder)
        def handler = Spy(new K8sTaskHandler(client: client, task:task))

        def POD_NAME = 'new-pod-id'
        def REQUEST =  [foo: 'bar']
        def RESPONSE = new K8sResponseJson([metadata: [name:POD_NAME]])
        def YAML = Paths.get('file.yaml')
        when:
        handler.submit()
        then:
        1 * handler.createBashWrapper(task) >> builder
        1 * builder.build() >> null
        1 * handler.yamlDebugPath() >> YAML
        1 * handler.newSubmitRequest(task) >> REQUEST
        1 * client.podCreate(REQUEST,YAML) >> RESPONSE
        handler.podName == POD_NAME
        handler.status == TaskStatus.SUBMITTED

        when:
        handler.submit()
        then:
        1 * handler.createBashWrapper(task) >> builder
        1 * builder.build() >> null
        1 * handler.yamlDebugPath() >> YAML
        1 * handler.newSubmitRequest(task) >> REQUEST
        1 * client.podCreate(REQUEST,YAML) >> new K8sResponseJson([missing: 'meta'])
        then:
        thrown(K8sResponseException)
    }

    def 'should submit a job' () {
        given:
        def WORK_DIR = Paths.get('/some/work/dir')
        def task = Mock(TaskRun)
        def client = Mock(K8sClient)
        def builder = Mock(K8sWrapperBuilder)
        def config = Mock(TaskConfig)
        def executor = Mock(K8sExecutor)
        def handler = Spy(new K8sTaskHandler(builder: builder, client: client, executor: executor))
        def podOptions = Mock(PodOptions)
        and:
        podOptions.automountServiceAccountToken >> true
        Map result

        when:
        result = handler.newSubmitRequest(task)
        then:
        1 * client.getConfig() >> new ClientConfig()
        _ * handler.fusionEnabled() >> false
        1 * handler.fixOwnership() >> false
        1 * handler.useJobResource() >> true
        1 * handler.entrypointOverride() >> true
        1 * handler.getSyntheticPodName(task) >> 'nf-123'
        1 * handler.getLabels(task) >> [:]
        1 * handler.getAnnotations() >> [:]
        1 * handler.getContainerMounts() >> []
        1 * handler.getPodOptions() >> podOptions
        1 * task.getContainer() >> 'debian:latest'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> config

        result == [
            apiVersion: 'batch/v1', 
            kind: 'Job', 
            metadata:[name: 'nf-123', namespace: 'default'], 
            spec:[
              backoffLimit: 0,
              template: [
                  spec: [
                     restartPolicy: 'Never',
                     containers: [
                       [
                           name: 'nf-123',
                           image: 'debian:latest',
                           command: ['/bin/bash', '-ue','/some/work/dir/.command.run']
                       ]
                     ]
                  ]
              ]
            ]
        ]
    }

    def 'should check if running'  () {
        given:
        def POD_NAME = 'pod-xyz'
        def client = Mock(K8sClient)
        def handler = Spy(new K8sTaskHandler(client: client, podName: POD_NAME, status: TaskStatus.SUBMITTED))

        when:
        def result = handler.checkIfRunning()
        then:
        1 * handler.getState() >> [:]
        result == false

        when:
        result = handler.checkIfRunning()
        then:
        1 * handler.getState() >> null
        result == false

        when:
        result = handler.checkIfRunning()
        then:
        1 * handler.getState() >> [running:["startedAt": "2018-01-13T10:19:16Z"]]
        result == true
    }

    def 'should check if completed' () {
        given:
        def ERR_FILE = Paths.get('err.file')
        def OUT_FILE = Paths.get('out.filex')
        def POD_NAME = 'pod-xyz'
        def EXIT_STATUS = 111
        def task = new TaskRun()
        def client = Mock(K8sClient)
        def termState = [ reason: "Completed",
                          startedAt: "2018-01-13T10:09:36Z",
                          finishedAt: "2018-01-13T10:19:36Z",
                          exitCode: 0 ]
        def fullState = [terminated: termState]
        and:
        def handler = Spy(new K8sTaskHandler(task: task, client:client, podName: POD_NAME, outputFile: OUT_FILE, errorFile: ERR_FILE))

        when:
        def result = handler.checkIfCompleted()
        then:
        1 * handler.getState() >> [:]
        handler.status != TaskStatus.COMPLETED
        result == false

        when:
        result = handler.checkIfCompleted()
        then:
        1 * handler.getState() >> null
        handler.status != TaskStatus.COMPLETED
        result == false

        when:
        result = handler.checkIfCompleted()
        then:
        1 * handler.getState() >> fullState
        1 * handler.updateTimestamps(termState)
        1 * handler.readExitFile() >> EXIT_STATUS
        1 * handler.deletePodIfSuccessful(task) >> null
        1 * handler.savePodLogOnError(task) >> null
        handler.task.exitStatus == EXIT_STATUS
        handler.task.@stdout == OUT_FILE
        handler.task.@stderr == ERR_FILE
        handler.status == TaskStatus.COMPLETED
        handler.startTimeMillis == 1515838176000
        handler.completeTimeMillis == 1515838776000
        result == true

    }

    def 'should kill a pod' () {
        given:
        def POD_NAME = 'pod-xyz'
        def client = Mock(K8sClient)
        def handler = Spy(new K8sTaskHandler(client:client, podName: POD_NAME))

        when:
        handler.kill()
        then:
        1 * handler.cleanupDisabled() >> false
        1 * client.podDelete(POD_NAME) >> null

        when:
        handler.kill()
        then:
        1 * handler.cleanupDisabled() >> true
        0 * client.podDelete(POD_NAME) >> null
    }

    def 'should check task cached state' () {
        given:
        def POD_NAME = 'pod-xyz'
        def client = Mock(K8sClient)
        def handler = Spy(new K8sTaskHandler(client:client, podName: POD_NAME))
        and:
        Map STATE1 = [status:'pending']
        Map STATE2 = [status:'running']
        Map STATE3 = [status:'complete']
        Map state

        // first time `client.podState` is invoked
        when:
        state = handler.getState()
        then:
        1 * client.podState(POD_NAME) >> STATE1
        state == STATE1

        // second time `client.podState` NOT invoked
        // the cached status is returned
        when:
        state = handler.getState()
        then:
        0 * client.podState(POD_NAME)
        state == STATE1

        // after more than a second `client.podState` is invoked
        // an empty value is returned, therefore the previous status is returned
        when:
        sleep 1_500
        state = handler.getState()
        then:
        1 * client.podState(POD_NAME) >> [:]
        state == STATE1

        // still an empty status
        // the previous status is returned
        when:
        state = handler.getState()
        then:
        1 * client.podState(POD_NAME) >> [:]
        state == STATE1

        // now, the a valid status is returned
        // the status is cached
        when:
        state = handler.getState()
        then:
        1 * client.podState(POD_NAME) >> STATE2
        state == STATE2

        // following invocation is cached
        // the previous status is returned
        when:
        state = handler.getState()
        then:
        0 * client.podState(POD_NAME)
        state == STATE2

        // after a second, the a new invocation is executed
        // the new status is returned
        when:
        sleep 1_500
        state = handler.getState()
        then:
        1 * client.podState(POD_NAME) >> STATE3
        state == STATE3
    }

    def 'should return nodeTermination state' () {
        given:
        def POD_NAME = 'pod-xyz'
        def client = Mock(K8sClient)
        def handler = Spy(new K8sTaskHandler(client:client, podName: POD_NAME))
        
        when:
        def state = handler.getState()
        then:
        1 * client.podState(POD_NAME) >> { throw new NodeTerminationException("Node shutdown happened") }
        then:
        state.terminated.startedAt
        state.terminated.finishedAt
        and:
        state.nodeTermination instanceof NodeTerminationException
        state.nodeTermination.message == "Node shutdown happened"
    }

    def 'should return container mounts' () {

        given:
        def wrapper = Mock(K8sWrapperBuilder)
        def k8sConfig = Mock(K8sConfig)
        and:
        def handler = Spy(new K8sTaskHandler(builder: wrapper))
        handler.getK8sConfig() >> k8sConfig

        when:
        def mounts = handler.getContainerMounts()
        then:
        1 * k8sConfig.getAutoMountHostPaths() >> false
        mounts == []

        when:
        mounts = handler.getContainerMounts()
        then:
        1 * k8sConfig.getAutoMountHostPaths() >> true
        1 * wrapper.getInputFiles() >> ['foo': Paths.get('/base_path/foo.txt'), 'bar': Paths.get('/base_path/bar.txt')]
        1 * wrapper.getBinDirs() >> [ Paths.get('/user/bin') ]
        1 * wrapper.getWorkDir() >> Paths.get('/work/dir')
        mounts == ['/base_path', '/user/bin', '/work/dir']

    }

    def 'should return labels map' () {
        given:
        def uuid = UUID.randomUUID()
        def task = Mock(TaskRun)
        def exec = Mock(K8sExecutor)
        def proc = Mock(TaskProcessor)
        def sess = Mock(Session)
        def handler = Spy(new K8sTaskHandler(executor: exec))

        when:
        def labels = handler.getLabels(task)
        then:
        handler.getRunName() >> 'pedantic-joe'
        task.getName() >> 'hello-world-1'
        task.getProcessor() >> proc
        task.getConfig() >> Mock(TaskConfig)
        proc.getName() >> 'hello-proc'
        exec.getSession() >> sess
        sess.getUniqueId() >> uuid
        exec.getK8sConfig() >> [pod: [
                [label: 'foo', value: 'bar'],
                [label: 'app', value: 'nextflow'],
                [label: 'x', value: 'hello_world']
        ]]

        labels.app == 'nextflow'
        labels.foo == 'bar'
        labels.x == 'hello_world'
        and:
        labels.'nextflow.io/app' == 'nextflow'
        labels.'nextflow.io/processName' == 'hello-proc'
        labels.'nextflow.io/taskName' ==  'hello-world-1'
        labels.'nextflow.io/sessionId' instanceof String
        labels.'nextflow.io/sessionId' == "uuid-${uuid.toString()}".toString()
        and:
        !labels.containsKey('nextflow.io/queue')
    }

    def 'should return process queue as a label'() {
        given:
        def uuid = UUID.randomUUID()
        def task = Mock(TaskRun)
        def exec = Mock(K8sExecutor)
        def proc = Mock(TaskProcessor)
        def sess = Mock(Session)
        def handler = Spy(new K8sTaskHandler(executor: exec))

        when:
        def labels = handler.getLabels(task)
        then:
        handler.getRunName() >> 'pedantic-joe'
        task.getName() >> 'hello-world-1'
        task.getProcessor() >> proc
        task.getConfig() >> new TaskConfig(queue: 'him-mem-queue')
        proc.getName() >> 'hello-proc'
        exec.getSession() >> sess
        sess.getUniqueId() >> uuid
        exec.getK8sConfig() >> [:]
        and:
        labels.'nextflow.io/queue' == 'him-mem-queue'
        and:
        labels.'nextflow.io/app' == 'nextflow'
        labels.'nextflow.io/processName' == 'hello-proc'
        labels.'nextflow.io/taskName' ==  'hello-world-1'
        labels.'nextflow.io/sessionId' instanceof String
        labels.'nextflow.io/sessionId' == "uuid-${uuid.toString()}".toString()
    }

    def 'should delete pod if complete' () {

        given:
        def POD_NAME = 'the-pod-name'
        def executor = Mock(K8sExecutor)
        def client = Mock(K8sClient)
        def handler = Spy(new K8sTaskHandler(podName: POD_NAME, executor:executor, client:client))
        and:
        def TASK_OK = Mock(TaskRun); TASK_OK.isSuccess() >> true
        def TASK_FAIL = Mock(TaskRun); TASK_FAIL.isSuccess() >> false

        when:
        handler.deletePodIfSuccessful(TASK_OK)
        then:
        1 * executor.getK8sConfig() >> new K8sConfig()
        1 * client.podDelete(POD_NAME) >> null

        when:
        handler.deletePodIfSuccessful(TASK_OK)
        then:
        1 * executor.getK8sConfig() >> new K8sConfig(cleanup: true)
        1 * client.podDelete(POD_NAME) >> null

        when:
        handler.deletePodIfSuccessful(TASK_FAIL)
        then:
        1 * executor.getK8sConfig() >> new K8sConfig(cleanup: false)
        0 * client.podDelete(POD_NAME) >> null

    }

    def 'should save pod log' () {

        given:
        def folder = Files.createTempDirectory('test')
        def POD_NAME = 'the-pod-name'
        def POD_MESSAGE = 'Hello world!'
        def POD_LOG = new ByteArrayInputStream(new String(POD_MESSAGE).bytes)
        def session = Mock(Session)
        def task = Mock(TaskRun)
        def executor = Mock(K8sExecutor)
        def client = Mock(K8sClient)
        and:
        def handler = Spy(new K8sTaskHandler(executor: executor, client: client, podName: POD_NAME))

        when:
        handler.savePodLogOnError(task)
        then:
        task.isSuccess() >> true
        0 * client.podLog(_)

        when:
        handler.savePodLogOnError(task)
        then:
        task.isSuccess() >> false
        task.getWorkDir() >> folder
        executor.getSession() >> session
        session.isTerminated() >> false
        session.isCancelled() >> false
        session.isAborted() >> false
        1 * client.podLog(POD_NAME) >> POD_LOG

        folder.resolve( TaskRun.CMD_LOG ).text == POD_MESSAGE
        cleanup:
        folder?.deleteDir()

    }

    def 'should merge pod options' () {

        given:
        PodOptions opts

        def taskConfig = Mock(TaskConfig)
        def task = Mock(TaskRun)
        task.getConfig() >> taskConfig

        def k8sConfig = Mock(K8sConfig)
        def handler = Spy(K8sTaskHandler)
        handler.getK8sConfig() >> k8sConfig
        handler.task = task

        when:
        opts = handler.getPodOptions()
        then:
        1 * taskConfig.getPodOptions() >> new PodOptions()
        1 * k8sConfig.getPodOptions() >> new PodOptions()
        opts == new PodOptions()

        when:
        opts = handler.getPodOptions()
        then:
        1 * taskConfig.getPodOptions() >> new PodOptions([[env:'HELLO', value:'WORLD']])
        1 * k8sConfig.getPodOptions() >> new PodOptions()
        opts == new PodOptions([[env:'HELLO', value:'WORLD']])


        when:
        opts = handler.getPodOptions()
        then:
        1 * taskConfig.getPodOptions() >> new PodOptions([[env:'HELLO', value:'WORLD']])
        1 * k8sConfig.getPodOptions() >> new PodOptions([ [env:'BRAVO', value:'HOTEL'] ])
        opts == new PodOptions([[env:'HELLO', value:'WORLD'], [env:'BRAVO', value:'HOTEL']])

        when:
        opts = handler.getPodOptions()
        then:
        1 * k8sConfig.getPodOptions() >> new PodOptions([[env: 'NXF_FUSION_BUCKETS', value: 's3://nextflow-ci'], [privileged: true]])
        and:
        1 * taskConfig.getPodOptions() >> new PodOptions([:])
        and: 
        opts == new PodOptions([[env: 'NXF_FUSION_BUCKETS', value: 's3://nextflow-ci'], [privileged: true]])
    }

    def 'should update startTimeMillis and completeTimeMillis with terminated state' () {

        given:
        def handler = Spy(K8sTaskHandler)
        def termState = [ startedAt: "2018-01-13T10:09:36Z",
                          finishedAt: "2018-01-13T10:19:36Z" ]

        when:
        handler.updateTimestamps(termState)
        then:
        handler.startTimeMillis == 1515838176000
        handler.completeTimeMillis == 1515838776000
    }

    def 'should update timestamps with current time with missing or malformed time' () {

        given:
        def handler = Spy(K8sTaskHandler)
        def malformedTime = [ startedAt: "2018-01-13 10:09:36",
                              finishedAt: "2018-01-13T10:19:36Z" ]

        def garbage = [ what: "nope" ]

        when:
        handler.updateTimestamps(malformedTime)
        then:
        handler.startTimeMillis > 0 // confirms that timestamps have been updated
        handler.startTimeMillis <= handler.completeTimeMillis // confirms that order is sane

        when:
        handler.updateTimestamps(garbage)
        then:
        handler.startTimeMillis > 0
        handler.startTimeMillis <= handler.completeTimeMillis
    }

    def 'should not update timestamps with malformed time and when startTimeMillis already set' () {

        given:
        def handler = Spy(K8sTaskHandler)
        handler.startTimeMillis = 10
        handler.completeTimeMillis = 20
        def malformedTime = [ startedAt: "2018-01-13 10:09:36",
                              finishedAt: "2018-01-13T10:19:36Z" ]

        when:
        handler.updateTimestamps(malformedTime)
        then:
        handler.startTimeMillis == 10
        handler.completeTimeMillis == 20
    }

    def 'should create a fusion pod' () {
        given:
        def WORK_DIR = XPath.get('http://some/work/dir')
        def config = Mock(TaskConfig)
        def task = Mock(TaskRun)
        def client = Mock(K8sClient)
        def builder = Mock(K8sWrapperBuilder)
        def handler = Spy(new K8sTaskHandler(builder:builder, client: client))
        Map result

        when:
        result = handler.newSubmitRequest(task)
        then:
        handler.getTask() >> task
        handler.fusionEnabled() >> true
        handler.fusionLauncher() >> Mock(FusionScriptLauncher) {
            fusionEnv() >> [NXF_FUSION_BUCKETS: 'this,that']
        }
        and:
        task.getContainer() >> 'debian:latest'
        task.getWorkDir() >> WORK_DIR
        task.getConfig() >> config
        and:
        1 * handler.fixOwnership() >> false
        1 * handler.entrypointOverride() >> false
        1 * handler.getPodOptions() >> new PodOptions()
        1 * handler.getSyntheticPodName(task) >> 'nf-123'
        1 * handler.getLabels(task) >> [:]
        1 * handler.getAnnotations() >> [:]
        1 * handler.getContainerMounts() >> []
        and:
        1 * config.getCpus() >> 0
        1 * config.getMemory() >> null
        1 * client.getConfig() >> new ClientConfig()
        result == [ apiVersion: 'v1',
                    kind: 'Pod',
                    metadata: [
                            name:'nf-123',
                            namespace:'default'
                    ],
                    spec: [
                            restartPolicy:'Never',
                            containers:[
                                    [name:'nf-123',
                                     image:'debian:latest',
                                     args:['bash', '-o', 'pipefail', '-c', 'trap "{ ret=$?; cp .command.log null||true; exit $ret; }" EXIT; bash null 2>&1 | tee .command.log'],
                                     securityContext:[privileged:true],
                                     env:[[name:'NXF_FUSION_BUCKETS', value:'this,that']]]
                            ]
                    ]
        ]
     }

    def 'get fusion submit command' () {
        given:
        def handler = Spy(K8sTaskHandler) {
            fusionEnabled() >> true
            fusionLauncher() >> new FusionScriptLauncher(scheme: 'http')
            getTask() >> Mock(TaskRun) {
                getWorkDir() >> XPath.get('http://foo/work/dir')
            }
        }

        when:
        def result =  handler.getSubmitCommand(Mock(TaskRun))
        then:
        result.join(' ') == 'bash -o pipefail -c trap "{ ret=$?; cp .command.log /fusion/http/foo/work/dir/.command.log||true; exit $ret; }" EXIT; bash /fusion/http/foo/work/dir/.command.run 2>&1 | tee .command.log'
    }
}
