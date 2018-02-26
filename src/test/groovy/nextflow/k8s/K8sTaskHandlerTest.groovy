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
import java.nio.file.Paths

import nextflow.Session
import nextflow.k8s.client.ClientConfig
import nextflow.k8s.client.K8sClient
import nextflow.k8s.client.K8sResponseException
import nextflow.k8s.client.K8sResponseJson
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sTaskHandlerTest extends Specification {


    def 'should return a new pod request with no storage' () {
        given:
        def WORK_DIR = Paths.get('/some/work/dir')
        def config = Mock(TaskConfig)
        def task = Mock(TaskRun)
        def client = Mock(K8sClient)
        def builder = Mock(K8sWrapperBuilder)
        def handler = Spy(K8sTaskHandler)
        handler.builder = builder
        handler.client = client
        Map result

        when:
        result = handler.newSubmitRequest(task)
        then:
        1 * handler.getSyntheticPodName(task) >> 'nf-123'
        1 * handler.getAutoMountHostPaths() >> false
        1 * handler.getVolumeClaims() >> [:]
        1 * handler.getLabels(task) >> [foo: 'bar', hello: 'world']
        0 * handler.getContainerMounts() >> null
        1 * task.getContainer() >> 'debian:latest'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> config
        1 * config.getCpus() >> 1
        1 * config.getMemory() >> null
        1 * client.getConfig() >> new ClientConfig()
        result == [ apiVersion: 'v1',
                    kind: 'Pod',
                    metadata: [
                            name:'nf-123',
                            namespace:'default',
                            labels:[ foo:'bar', hello: 'world'] ],
                    spec: [
                            restartPolicy:'Never',
                            containers:[
                                    [name:'nf-123',
                                     image:'debian:latest',
                                     command:['/bin/bash', '-ue','.command.run'],
                                     workingDir:'/some/work/dir']
                            ]
                    ]
                ]

        when:
        result = handler.newSubmitRequest(task)
        then:
        1 * handler.getSyntheticPodName(task) >> 'nf-foo'
        1 * handler.getAutoMountHostPaths() >> false
        1 * handler.getVolumeClaims() >> null
        1 * handler.getLabels(task) >> [sessionId:'xxx']
        0 * handler.getContainerMounts() >> null
        1 * builder.fixOwnership() >> true
        1 * handler.getOwner() >> '501:502'
        1 * task.getContainer() >> 'debian:latest'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> config
        1 * config.getCpus() >> 1
        1 * config.getMemory() >> null
        1 * client.getConfig() >> new ClientConfig()
        result == [ apiVersion: 'v1',
                    kind: 'Pod',
                    metadata: [name:'nf-foo', namespace:'default', labels: [sessionId: 'xxx']],
                    spec: [
                            restartPolicy:'Never',
                            containers:[
                                    [name:'nf-foo',
                                     image:'debian:latest',
                                     command:['/bin/bash', '-ue','.command.run'],
                                     workingDir:'/some/work/dir',
                                     env: [ [name:'NXF_OWNER', value:'501:502'] ]
                                    ]
                            ]
                    ]
        ]


        when:
        result = handler.newSubmitRequest(task)
        then:
        1 * handler.getSyntheticPodName(task) >> 'nf-abc'
        1 * handler.getLabels(task) >> null
        1 * handler.getAutoMountHostPaths() >> false
        1 * handler.getVolumeClaims() >> null
        0 * handler.getContainerMounts() >> null
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
                                     command:['/bin/bash', '-ue', '.command.run'],
                                     workingDir:'/some/work/dir',
                                     resources:[ limits:[cpu:4, memory:'16384Mi'] ]
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
        def handler = Spy(K8sTaskHandler)
        def config = Mock(ClientConfig)
        handler.builder = builder
        handler.client = client
        Map result

        when:
        result = handler.newSubmitRequest(task)
        then:
        1 * handler.getSyntheticPodName(task) >> 'nf-123'
        1 * handler.getAutoMountHostPaths() >> false
        1 * handler.getVolumeClaims() >> [:]
        1 * handler.getLabels(task) >> [:]
        0 * handler.getContainerMounts() >> null
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
                                     command:['/bin/bash', '-ue','.command.run'],
                                     workingDir:'/some/work/dir']
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
        def handler = Spy(K8sTaskHandler)
        handler.builder = builder
        handler.client = client
        Map result

        def CLAIMS = [ 'first': [mountPath: '/work'], second: [mountPath: '/data'] ]

        when:
        result = handler.newSubmitRequest(task)
        then:
        1 * handler.getSyntheticPodName(task) >> 'nf-123'
        1 * handler.getVolumeClaims() >> CLAIMS
        1 * handler.getAutoMountHostPaths() >> false
        0 * handler.getContainerMounts() >> null
        1 * handler.getLabels(task) >> null
        1 * task.getContainer() >> 'debian:latest'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> config
        1 * config.getCpus() >> 1
        1 * config.getMemory() >> null
        1 * client.getConfig() >> new ClientConfig()

        result == [ apiVersion: 'v1',
                    kind: 'Pod',
                    metadata: [name:'nf-123', namespace:'default' ],
                    spec: [
                            restartPolicy:'Never',
                            containers:[
                                    [name: 'nf-123',
                                     image: 'debian:latest',
                                     command: ['/bin/bash', '-ue', '.command.run'],
                                     workingDir: '/some/work/dir',
                                     volumeMounts: [
                                             [name:'vol-1', mountPath:'/work'],
                                             [name:'vol-2', mountPath:'/data']
                                     ] ]
                            ],
                            volumes: [   [name:'vol-1', persistentVolumeClaim:[claimName: 'first']],
                                         [name:'vol-2', persistentVolumeClaim:[claimName: 'second']] ]
                    ]
        ]


        when:
        result = handler.newSubmitRequest(task)
        then:
        1 * handler.getSyntheticPodName(task) >> 'nf-123'
        1 * handler.getVolumeClaims() >> null
        1 * handler.getAutoMountHostPaths() >> true
        1 * handler.getContainerMounts() >> ['/tmp', '/data']
        1 * handler.getLabels(task) >> null
        1 * task.getContainer() >> 'debian:latest'
        1 * task.getWorkDir() >> WORK_DIR
        1 * task.getConfig() >> config
        1 * config.getCpus() >> 1
        1 * config.getMemory() >> null
        1 * client.getConfig() >> new ClientConfig()

        result == [ apiVersion: 'v1',
                    kind: 'Pod',
                    metadata: [name:'nf-123', namespace:'default' ],
                    spec: [
                            restartPolicy:'Never',
                            containers:[
                                    [name: 'nf-123',
                                     image: 'debian:latest',
                                     command: ['/bin/bash', '-ue', '.command.run'],
                                     workingDir: '/some/work/dir',
                                     volumeMounts: [
                                             [name:'vol-3', mountPath:'/tmp'],
                                             [name:'vol-4', mountPath: '/data']
                                     ] ]
                            ],
                            volumes: [  [name:'vol-3', hostPath:[path:'/tmp']],
                                        [name:'vol-4', hostPath:[path:'/data']]]
                    ]
        ]

    }

    def 'should submit a pod' () {

        given:
        def task = Mock(TaskRun)
        def client = Mock(K8sClient)
        def builder = Mock(K8sWrapperBuilder)
        def handler = Spy(K8sTaskHandler)
        handler.client = client
        handler.task = task

        final POD_NAME = 'new-pod-id'
        final REQUEST =  [foo: 'bar']
        final RESPONSE = new K8sResponseJson([metadata: [name:POD_NAME]])
        final YAML = Paths.get('file.yaml')
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



    def 'should check if running'  () {
        given:
        final POD_NAME = 'pod-xyz'
        def client = Mock(K8sClient)
        def handler = Spy(K8sTaskHandler)
        handler.client = client
        handler.podName = POD_NAME

        when:
        def result = handler.checkIfRunning()
        then:
        1 * handler.getState() >> [:]
        result == false

        when:
        result = handler.checkIfRunning()
        then:
        1 * handler.getState() >> [running:["startedAt": "2018-01-13T10:19:16Z"]]
        result == true
    }

    def 'should check if completed' () {
        given:
        final ERR_FILE = Paths.get('err.file')
        final OUT_FILE = Paths.get('out.filex')
        final POD_NAME = 'pod-xyz'
        final EXIT_STATUS = 111
        def task = new TaskRun()
        def client = Mock(K8sClient)
        def handler = Spy(K8sTaskHandler)
        handler.task = task
        handler.client = client
        handler.podName = POD_NAME
        handler.outputFile = OUT_FILE
        handler.errorFile = ERR_FILE

        when:
        def result = handler.checkIfCompleted()
        then:
        1 * handler.getState() >> [:]
        handler.status != TaskStatus.COMPLETED
        result == false

        when:
        result = handler.checkIfCompleted()
        then:
        1 * handler.getState() >> [terminated: ["reason": "Completed", "finishedAt": "2018-01-13T10:19:36Z", exitCode: 0]]
        1 * handler.readExitFile() >> EXIT_STATUS
        1 * handler.deletePodIfSuccessful(task) >> null
        handler.task.exitStatus == EXIT_STATUS
        handler.task.@stdout == OUT_FILE
        handler.task.@stderr == ERR_FILE
        handler.status == TaskStatus.COMPLETED
        result == true

    }

    def 'should kill a pod' () {
        given:
        final POD_NAME = 'pod-xyz'
        def client = Mock(K8sClient)
        def handler = Spy(K8sTaskHandler)
        handler.client = client
        handler.podName = POD_NAME

        when:
        handler.kill()
        then:
        1 * client.podDelete(POD_NAME) >> null
    }

    def 'should check task cached state' () {
        given:
        final POD_NAME = 'pod-xyz'
        def client = Mock(K8sClient)
        def handler = Spy(K8sTaskHandler)
        handler.client = client
        handler.podName = POD_NAME
        Map STATE1 = [foo:1, bar:2]
        Map STATE2 = [fizz:1, bazz: 2]
        Map state

        when:
        state = handler.getState()
        then:
        1 * client.podState(POD_NAME) >> STATE1
        state == STATE1

        when:
        state = handler.getState()
        then:
        0 * client.podState(POD_NAME) >> 0
        state == STATE1

        when:
        sleep 1_500
        state = handler.getState()
        then:
        1 * client.podState(POD_NAME) >> STATE2
        state == STATE2
    }

    def 'should return container mounts' () {

        given:
        def wrapper = Mock(K8sWrapperBuilder)
        def handler = new K8sTaskHandler(builder: wrapper)

        when:
        def mounts = handler.getContainerMounts()
        then:
        1 * wrapper.getResolvedInputs() >> ['foo': Paths.get('/base_path/foo.txt'), 'bar': Paths.get('/base_path/bar.txt')]
        1 * wrapper.getBinDir() >> Paths.get('/user/bin')
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
        def handler = Spy(K8sTaskHandler)
        handler.executor = exec

        when:
        def labels = handler.getLabels(task)
        then:
        handler.getRunName() >> 'pedantic-joe'
        task.getName() >> 'hello-world-1'
        task.getProcessor() >> proc 
        proc.getName() >> 'hello-proc'
        exec.getSession() >> sess
        sess.getUniqueId() >> uuid
        exec.getK8sConfig() >> [labels: [foo:'bar', app: 'snakemake', x:'hello world', y:null ]]

        labels.app == 'nextflow'
        labels.foo == 'bar'
        labels.x == 'hello_world'    
        labels.y == 'null'
        labels.processName == 'hello-proc'
        labels.taskName ==  'hello-world-1'
        labels.sessionId instanceof String 
        labels.sessionId == "uuid-${uuid.toString()}".toString()
    }

    @Unroll
    def 'should sanitize k8s label: #label' () {

        given:
        def handler = new K8sTaskHandler()

        expect:
        handler.sanitize0(label) == str

        where:
        label           | str
        null            | 'null'
        'hello'         | 'hello'
        'hello world'   | 'hello_world'
        'hello  world'  | 'hello_world'
        'hello.world'   | 'hello.world'
        'hello-world'   | 'hello-world'
        'hello_world'   | 'hello_world'
        'hello_world-'  | 'hello_world'
        'hello_world_'  | 'hello_world'
        'hello_world.'  | 'hello_world'
        'hello_123'     | 'hello_123'
        'HELLO 123'     | 'HELLO_123'
        '123hello'      | 'hello'
    }

    def 'should delete pod if complete' () {

        given:
        def POD_NAME = 'the-pod-name'
        def executor = Mock(K8sExecutor)
        def client = Mock(K8sClient)
        def handler = Spy(K8sTaskHandler)
        handler.podName = POD_NAME
        handler.executor = executor
        handler.client = client
        def TASK_OK = Mock(TaskRun); TASK_OK.isSuccess() >> true
        def TASK_FAIL = Mock(TaskRun); TASK_FAIL.isSuccess() >> false

        when:
        handler.deletePodIfSuccessful(TASK_OK)
        then:
        1 * executor.getK8sConfig() >> [cleanup: true]
        1 * client.podDelete(POD_NAME) >> null

        when:
        handler.deletePodIfSuccessful(TASK_OK)
        then:
        1 * executor.getK8sConfig() >> [cleanup: false]
        0 * client.podDelete(POD_NAME) >> null


        when:
        handler.deletePodIfSuccessful(TASK_FAIL)
        then:
        1 * executor.getK8sConfig() >> [cleanup: true]
        0 * client.podDelete(POD_NAME) >> null

    }
}
