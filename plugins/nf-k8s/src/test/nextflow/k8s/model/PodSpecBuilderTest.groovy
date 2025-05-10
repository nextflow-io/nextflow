/*
 * Copyright 2013-2024, Seqera Labs
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

import nextflow.executor.res.AcceleratorResource
import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PodSpecBuilderTest extends Specification {

    def setup() {
        PodSpecBuilder.VOLUMES.set(0)
    }


    def 'should create pod spec' () {

        when:
        def spec = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withWorkDir('/some/work/dir')
                .withCommand(['echo', 'hello'])
                .build()

        then:
        spec == [
                apiVersion: 'v1',
                kind: 'Pod',
                metadata: [name:'foo', namespace:'default'],
                spec: [
                        restartPolicy:'Never',
                        containers:[[
                                name:'foo',
                                image:'busybox',
                                command:['echo', 'hello'],
                                workingDir:'/some/work/dir'
                        ]]
                ]
        ]

    }

    def 'should create pod spec with args' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withArgs(['echo', 'hello'])
                .build()

        then:
        pod.spec.containers[0].args == ['echo', 'hello']

    }

    def 'should create pod spec with args string' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withArgs('echo foo')
                .build()

        then:
        pod.spec.containers[0].args == ['/bin/bash', '-c', 'echo foo']

    }

    def 'should create pod spec with privileged' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand('echo foo')
                .withPrivileged(true)
                .build()

        then:
        pod.spec.containers[0].securityContext == [privileged: true]

    }

    def 'should create pod with resources limits' () {
        when:
        def pod1 = new PodSpecBuilder()
            .withPodName('foo')
            .withImageName('busybox')
            .withCommand('echo foo')
            .withResourcesLimits('nextflow.io/fuse': 1)
            .build()

        then:
        pod1.spec.containers[0].resources == [limits:['nextflow.io/fuse':1]]


        when:
        def pod2 = new PodSpecBuilder()
            .withPodName('foo')
            .withImageName('busybox')
            .withCommand('echo foo')
            .withCpus(8)
            .withCpuLimits(true)
            .withMemory(MemoryUnit.of('10GB'))
            .withResourcesLimits('nextflow.io/fuse': 1)
            .build()

        then:
        pod2.spec.containers[0].resources == [
                requests: ['cpu':8, 'memory':'10240Mi'],
                limits: ['cpu':8, 'memory':'10240Mi', 'nextflow.io/fuse':1] ]
    }

    def 'should set namespace, labels and annotations' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['sh', '-c', 'echo hello'])
                .withNamespace('xyz')
                .withLabel('app','myApp')
                .withLabel('runName','something')
                .withLabel('version','3.6.1')
                .withAnnotation("anno1", "value1")
                .withAnnotations([anno2: "value2", anno3: "value3"])
                .build()

        then:
        pod.metadata.namespace == 'xyz'
        pod.metadata.labels == [
                app: 'myApp',
                runName: 'something',
                version: '3.6.1'
        ]
        pod.metadata.annotations == [
                anno1: "value1",
                anno2: "value2",
                anno3: "value3"
        ]
    }
    
    def 'should truncate labels longer than 63 chars' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['sh', '-c', 'echo hello'])
                .withLabel('app','myApp')
                .withLabel('runName','something')
                .withLabel('tag','somethingreallylonggggggggggggggggggggggggggggggggggggggggggendEXTRABIT')
                .withLabels([tag2: 'somethingreallylonggggggggggggggggggggggggggggggggggggggggggendEXTRABIT', tag3: 'somethingreallylonggggggggggggggggggggggggggggggggggggggggggendEXTRABIT'])
                .build()

        then:
        pod.metadata.labels == [
                app: 'myApp',
                runName: 'something',
                tag: 'somethingreallylonggggggggggggggggggggggggggggggggggggggggggend',
                tag2: 'somethingreallylonggggggggggggggggggggggggggggggggggggggggggend',
                tag3: 'somethingreallylonggggggggggggggggggggggggggggggggggggggggggend'
        ]
    }


    def 'should set resources and env' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand('echo hello')
                .withEnv(PodEnv.value('ALPHA','hello'))
                .withEnv(PodEnv.value('DELTA', 'world'))
                .withCpus(8)
                .withAccelerator( new AcceleratorResource(request: 5, limit:10, type: 'foo.org') )
                .withMemory('100Gi')
                .withDisk('10Gi')
                .build()

        then:
        pod.spec.containers[0].env == [
                [name:'ALPHA', value:'hello'],
                [name:'DELTA', value:'world']
        ]
        pod.spec.containers[0].resources == [
                requests: ['foo.org/gpu':5, cpu:8, memory:'100Gi', 'ephemeral-storage':'10Gi'],
                limits: ['foo.org/gpu':10, memory:'100Gi', 'ephemeral-storage':'10Gi']
        ]
    }

    def 'should get storage spec for volume claims' () {

        when:
        def pod = new PodSpecBuilder()
                    .withPodName('foo')
                    .withImageName('busybox')
                    .withCommand(['echo'])
                    .withVolumeClaim(new PodVolumeClaim('first','/work'))
                    .withVolumeClaim(new PodVolumeClaim('second', '/data', '/foo'))
                    .withVolumeClaim(new PodVolumeClaim('third', '/things', null, true))
                    .build()
        then:
        pod.spec.containers[0].volumeMounts == [
                [name:'vol-1', mountPath:'/work'],
                [name:'vol-2', mountPath:'/data', subPath: '/foo'],
                [name:'vol-3', mountPath:'/things', readOnly: true]
        ]
        pod.spec.volumes == [
                [name:'vol-1', persistentVolumeClaim:[claimName:'first']],
                [name:'vol-2', persistentVolumeClaim:[claimName:'second']],
                [name:'vol-3', persistentVolumeClaim:[claimName:'third']]
        ]

    }

    def 'should only define one volume per persistentVolumeClaim' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo'])
                .withVolumeClaim(new PodVolumeClaim('first','/work'))
                .withVolumeClaim(new PodVolumeClaim('first','/work2', '/bar'))
                .withVolumeClaim(new PodVolumeClaim('second', '/data', '/foo'))
                .withVolumeClaim(new PodVolumeClaim('second', '/data2', '/fooz'))
                .build()
        then:
        pod.spec.containers[0].volumeMounts == [
                [name:'vol-1', mountPath:'/work'],
                [name:'vol-1', mountPath:'/work2', subPath: '/bar'],
                [name:'vol-2', mountPath:'/data', subPath: '/foo'],
                [name:'vol-2', mountPath:'/data2', subPath: '/fooz']
        ]
        pod.spec.volumes == [
                [name:'vol-1', persistentVolumeClaim:[claimName:'first']],
                [name:'vol-2', persistentVolumeClaim:[claimName:'second']]
        ]

    }

    def 'should get config map mounts' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo'])
                .withConfigMap(new PodMountConfig(config: 'cfg1', mountPath: '/etc/config'))
                .withConfigMap(new PodMountConfig(config: 'data2', mountPath: '/data/path'))
                .build()
        then:
        pod.spec.containers[0].volumeMounts == [
                [name:'vol-1', mountPath:'/etc/config'],
                [name:'vol-2', mountPath:'/data/path']
        ]
        pod.spec.volumes == [
                [name:'vol-1', configMap:[name:'cfg1']],
                [name:'vol-2', configMap:[name:'data2']]
        ]

    }

    def 'should get csi ephemeral mounts' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo'])
                .withCsiEphemeral(new PodMountCsiEphemeral(csi: [driver: 'inline.storage.kubernetes.io', readOnly: true], mountPath: '/data'))
                .build()
        then:
        pod.spec.containers[0].volumeMounts == [
                [name: 'vol-1', mountPath: '/data', readOnly: true]
        ]
        pod.spec.volumes == [
                [name: 'vol-1', csi: [driver: 'inline.storage.kubernetes.io', readOnly: true]]
        ]
    }

    def 'should get empty dir mounts' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo'])
                .withEmptyDir(new PodMountEmptyDir(mountPath: '/scratch1', emptyDir: [medium: 'Disk']))
                .withEmptyDir(new PodMountEmptyDir(mountPath: '/scratch2', emptyDir: [medium: 'Memory']))
                .build()
        then:
        pod.spec.containers[0].volumeMounts == [
                [name: 'vol-1', mountPath: '/scratch1'],
                [name: 'vol-2', mountPath: '/scratch2']
        ]
        pod.spec.volumes == [
                [name: 'vol-1', emptyDir: [medium: 'Disk']],
                [name: 'vol-2', emptyDir: [medium: 'Memory']]
        ]
    }

    def 'should consume env secrets' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo'])
                .withEnv( PodEnv.value('FOO','abc'))
                .withEnv( PodEnv.secret('VAR_X', 'delta/bar'))
                .withEnv( PodEnv.secret('VAR_Y', 'gamma'))
                .build()

        then:
        pod.spec.containers[0].env == [
                [name: 'FOO', value: 'abc'],
                [name: 'VAR_X', valueFrom: [secretKeyRef: [name:'delta', key:'bar']]],
                [name: 'VAR_Y', valueFrom: [secretKeyRef: [name:'gamma', key:'VAR_Y']]]
        ]
    }

    def 'should consume env configMap' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo'])
                .withEnv( PodEnv.value('FOO','abc'))
                .withEnv( PodEnv.config('VAR_X', 'data'))
                .withEnv( PodEnv.config('VAR_Y', 'omega/bar-2'))
                .build()

        then:
        pod.spec.containers[0].env == [
                [name: 'FOO', value: 'abc'],
                [name: 'VAR_X', valueFrom: [configMapKeyRef: [name:'data', key:'VAR_X']]],
                [name: 'VAR_Y', valueFrom: [configMapKeyRef: [name:'omega', key:'bar-2']]]
        ]
    }

    def 'should consume file secrets' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo'])
                .withSecret(new PodMountSecret(secret: 'alpha', mountPath: '/this/and/that'))
                .withSecret(new PodMountSecret(secret: 'delta/foo', mountPath: '/etc/mnt/bar.txt'))
                .build()

        then:
        pod.spec.containers[0].volumeMounts == [
                [name:'vol-1', mountPath:'/this/and/that'],
                [name:'vol-2', mountPath:'/etc/mnt']
        ]
        pod.spec.volumes == [
                [name:'vol-1', secret:[secretName: 'alpha']],
                [name:'vol-2', secret:[
                        secretName: 'delta',
                        items: [
                                [ key: 'foo', path:'bar.txt' ]
                        ]
                ]]
        ]
    }

    def 'should get host path mounts' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo'])
                .withHostMount('/tmp','/scratch')
                .withHostMount('/host/data','/mnt/container')
                .build()

        then:
        pod.spec.containers[0].volumeMounts == [
                [name:'vol-1', mountPath:'/scratch'],
                [name:'vol-2', mountPath:'/mnt/container']
        ]
        pod.spec.volumes == [
                [name:'vol-1', hostPath: [path:'/tmp']],
                [name:'vol-2', hostPath: [path:'/host/data']]
        ]

    }


    def 'should return secret file volume and mounts' () {

        given:
        List mounts
        List volumes
        def builder = new PodSpecBuilder()

        when:
        def secret1 = new PodMountSecret(secret:'foo', mountPath: '/etc/conf')
        builder.secretToSpec( 'vol1', secret1, mounts=[], volumes=[] )

        then:
        mounts == [
                [ name: 'vol1', mountPath: '/etc/conf']
        ]

        volumes == [
                [ name: 'vol1', secret: [secretName: 'foo']]
        ]


        when:
        def secret2 = new PodMountSecret(secret:'bar/hello.txt', mountPath: '/etc/conf/world.txt')
        builder.secretToSpec( 'vol2', secret2, mounts=[], volumes=[] )

        then:
        mounts == [
                [ name: 'vol2', mountPath: '/etc/conf']
        ]

        volumes == [
                [ name: 'vol2', secret: [
                        secretName: 'bar',
                        items: [ [key: 'hello.txt', path:'world.txt'] ]
                ]]
        ]

    }

    def 'should return configmap file volume and mounts' () {

        given:
        List mounts
        List volumes
        def builder = new PodSpecBuilder()

        when:
        def config1 = new PodMountConfig(config:'foo', mountPath: '/etc/conf')
        builder.configMapToSpec( 'vol1', config1, mounts=[], volumes=[] )

        then:
        mounts == [
                [ name: 'vol1', mountPath: '/etc/conf']
        ]

        volumes == [
                [ name: 'vol1', configMap: [name: 'foo']]
        ]


        when:
        def config2 = new PodMountConfig(config:'bar/hello.txt', mountPath: '/etc/conf/world.txt')
        builder.configMapToSpec( 'vol2', config2, mounts=[], volumes=[] )

        then:
        mounts == [
                [ name: 'vol2', mountPath: '/etc/conf']
        ]

        volumes == [
                [ name: 'vol2', configMap: [
                        name: 'bar',
                        items: [ [key: 'hello.txt', path:'world.txt'] ]
                ]]
        ]

    }


    def 'should create pod spec with pod options' () {

        given:
        def affinity = [
            nodeAffinity: [
                requiredDuringSchedulingIgnoredDuringExecution: [
                    nodeSelectorTerms: [
                        [key: 'foo', operator: 'In', values: ['bar', 'baz']]
                    ]
                ]
            ]
        ]
        def tolerations = [[
            key: 'example-key',
            operator: 'Exists',
            effect: 'NoSchedule'
        ]]
        def opts = Mock(PodOptions)
        and:
        def builder = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo'])
                .withLabel('runName', 'crazy_john')
                .withAnnotation('evict', 'false')

        when:
        def pod = builder.withPodOptions(opts).build()
        then:
        _ * opts.getAffinity() >> affinity
        _ * opts.getAnnotations() >> [OMEGA:'zzz', SIGMA:'www']
        _ * opts.getAutomountServiceAccountToken() >> false
        2 * opts.getEnvVars() >> [ PodEnv.value('HELLO','WORLD') ]
        2 * opts.getImagePullPolicy() >> 'always'
        2 * opts.getImagePullSecret() >> 'myPullSecret'
        _ * opts.getLabels() >> [ALPHA: 'xxx', GAMMA: 'yyy']
        2 * opts.getVolumeClaims() >> [ new PodVolumeClaim('pvc1', '/work') ]
        2 * opts.getMountConfigMaps() >> [ new PodMountConfig('data', '/home/user') ]
        2 * opts.getMountSecrets() >> [ new PodMountSecret('blah', '/etc/secret.txt') ]
        _ * opts.getNodeSelector() >> new PodNodeSelector(gpu:true, queue: 'fast')
        _ * opts.getPriorityClassName() >> 'high-priority'
        _ * opts.getSecurityContext() >> new PodSecurityContext(1000)
        _ * opts.getTolerations() >> tolerations
        and:
        pod.metadata == [
                name:'foo',
                namespace:'default',
                labels:[runName:'crazy_john', ALPHA:'xxx', GAMMA:'yyy'],
                annotations: [evict: 'false', OMEGA:'zzz', SIGMA:'www']
        ]
        and:
        pod.spec.affinity == affinity
        pod.spec.automountServiceAccountToken == false
        pod.spec.imagePullSecrets == [[ name: 'myPullSecret' ]]
        pod.spec.nodeSelector == [gpu: 'true', queue: 'fast']
        pod.spec.priorityClassName == 'high-priority'
        pod.spec.securityContext == [ runAsUser: 1000 ]
        pod.spec.tolerations == tolerations
        pod.spec.containers[0].imagePullPolicy == 'always'
        pod.spec.containers[0].env == [[name:'HELLO', value:'WORLD']]
        pod.spec.containers[0].volumeMounts == [
                [name:'vol-1', mountPath:'/work'],
                [name:'vol-2', mountPath:'/home/user'],
                [name:'vol-3', mountPath:'/etc/secret.txt']
        ]
        and:
        pod.spec.volumes == [
                [name:'vol-1', persistentVolumeClaim:[claimName:'pvc1']],
                [name:'vol-2', configMap:[name:'data']],
                [name:'vol-3', secret:[secretName:'blah']]
        ]

    }

    def 'should create pod spec with activeDeadlineSeconds' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo', 'hello'])
                .withActiveDeadline(100)
                .build()

        then:
        pod.spec.activeDeadlineSeconds == 100

    }

    def 'should create pod spec with schedulerName' () {

        when:
        def pod = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo', 'hello'])
                .withPodOptions(new PodOptions(schedulerName: 'my-scheduler'))
                .build()

        then:
        pod.spec.schedulerName == 'my-scheduler'

    }

    def 'should create image pull request map' () {
        given:
        def builder = new PodSpecBuilder(imagePullSecret: 'MySecret')
        when:
        def result = builder.createPullSecret()
        then:
        result.size() == 1 
        result.get(0).name == 'MySecret'
    }


    def 'should return the resources map' () {

        given:
        def builder = new PodSpecBuilder()

        when:
        def res = builder.addAcceleratorResources(new AcceleratorResource(request:2, limit: 5), null)
        then:
        res.requests == ['nvidia.com/gpu': 2]
        res.limits == ['nvidia.com/gpu': 5]

        when:
        res = builder.addAcceleratorResources(new AcceleratorResource(limit: 5, type:'foo'), null)
        then:
        res.requests == ['foo.com/gpu': 5]
        res.limits == ['foo.com/gpu': 5]

        when:
        res = builder.addAcceleratorResources(new AcceleratorResource(request: 5, type:'foo.org'), null)
        then:
        res.requests == ['foo.org/gpu': 5]
        res.limits == null

        when:
        res = builder.addAcceleratorResources(new AcceleratorResource(request: 5, type: 'foo.org'), [requests: [cpu: 2]])
        then:
        res.requests == [cpu: 2, 'foo.org/gpu': 5]
        res.limits == null

        when:
        res = builder.addAcceleratorResources(new AcceleratorResource(request: 5, limit: 10, type: 'foo.org'), [requests: [cpu: 2]])
        then:
        res.requests == [cpu: 2, 'foo.org/gpu': 5]
        res.limits == ['foo.org/gpu': 10]

        when:
        res = builder.addAcceleratorResources(new AcceleratorResource(request: 5, type:'example.com/fpga'), null)
        then:
        res.requests == ['example.com/fpga': 5]
        res.limits == null

        when:
        res = builder.addAcceleratorResources(new AcceleratorResource(request: 5, limit: 10, type: 'example.com/fpga'), [requests: [cpu: 2]])
        then:
        res.requests == [cpu: 2, 'example.com/fpga': 5]
        res.limits == ['example.com/fpga': 10]
    }

    def 'should add resources limits' () {
        given:
        def builder = new PodSpecBuilder()
        Map resources

        when:
        resources = builder.addResourcesLimits(['foo':1], null)
        then:
        resources == [limits:[foo:1]]

        when:
        resources = builder.addResourcesLimits(['foo':1], [requests: ['x':1], limits: ['y': 2]])
        then:
        resources == [requests:[x:1], limits:[y:2, foo:1]]
    }


    @Unroll
    def 'should sanitize k8s label value: #label' () {
        given:
        def builder = new PodSpecBuilder()

        expect:
        builder.sanitizeValue(label, PodSpecBuilder.MetaType.LABEL, PodSpecBuilder.SegmentType.VALUE) == str

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
        '123hello'      | '123hello'
        'x2345678901234567890123456789012345678901234567890123456789012345' | 'x23456789012345678901234567890123456789012345678901234567890123'
    }

    @Unroll
    def 'should sanitize k8s label key: #label_key' () {
        given:
        def builder = new PodSpecBuilder()

        expect:
        builder.sanitizeKey(label_key, PodSpecBuilder.MetaType.LABEL) == str

        where:
        label_key               | str
        'foo'                   | 'foo'
        'key 1'                 | 'key_1'
        'foo.bar/key 2'         | 'foo.bar/key_2'
        'foo.bar/'              | 'foo.bar'
        '/foo.bar'              | 'foo.bar'
        'x2345678901234567890123456789012345678901234567890123456789012345' | 'x23456789012345678901234567890123456789012345678901234567890123'
        'x23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345/key 2' | 'x234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123/key_2'
        'foo.bar/x2345678901234567890123456789012345678901234567890123456789012345' | 'foo.bar/x23456789012345678901234567890123456789012345678901234567890123'
    }

    @Unroll
    def 'should report error if sanitizing k8s label with more than one slash character: #label_key' () {
        given:
        def builder = new PodSpecBuilder()

        when:
        builder.sanitizeKey(label_key, PodSpecBuilder.MetaType.LABEL)

        then:
        def error = thrown(expectedException)
        error.message == expectedMessage

        where:
        label_key               | expectedException         | expectedMessage
        'foo.bar/key 2/key 3'   | IllegalArgumentException  | "Invalid key in pod label -- Key can only contain exactly one '/' character"
        'foo.bar/foo/bar/bar'   | IllegalArgumentException  | "Invalid key in pod label -- Key can only contain exactly one '/' character"
    }

    @Unroll
    def 'should sanitize k8s label map' () {
        given:
        def builder = new PodSpecBuilder()

        expect:
        builder.sanitize(KEY_VALUE, PodSpecBuilder.MetaType.LABEL) == EXPECTED

        where:
        KEY_VALUE                   | EXPECTED
        [foo:'bar']                 | [foo:'bar']
        ['key 1':'value 2']         | [key_1:'value_2']
        ['foo.bar/key 2':'value 3'] | ['foo.bar/key_2':'value_3']
    }

    @Unroll
    def 'should sanitize k8s annotation key' () {
        given:
        def builder = new PodSpecBuilder()

        expect:
        builder.sanitize(KEY_VALUE, PodSpecBuilder.MetaType.ANNOTATION) == EXPECTED

        where:
        KEY_VALUE                           | EXPECTED
        [foo:'bar']                         | [foo:'bar']
        ['key 1':'value 2']                 | [key_1:'value 2']
        ['foo.bar/key 2':'value 3']         | ['foo.bar/key_2':'value 3']
        ['x2345678901234567890123456789012345678901234567890123456789012345':'value 5'] | ['x23456789012345678901234567890123456789012345678901234567890123':'value 5']
        ['x23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345/key 4':'value 6'] | ['x234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123/key_4':'value 6']
        ['foo.bar/x2345678901234567890123456789012345678901234567890123456789012345':'value 7'] | ['foo.bar/x23456789012345678901234567890123456789012345678901234567890123':'value 7']
    }

    @Unroll
    def 'should report error if sanitizing k8s annotation key with more than one slash character: #annotation_key' () {
        given:
        def builder = new PodSpecBuilder()

        when:
        builder.sanitizeKey(annotation_key, PodSpecBuilder.MetaType.ANNOTATION)

        then:
        def error = thrown(expectedException)
        error.message == expectedMessage

        where:
        annotation_key          | expectedException         | expectedMessage
        'foo.bar/key 2/key 3'   | IllegalArgumentException  | "Invalid key in pod annotation -- Key can only contain exactly one '/' character"
        'foo.bar/foo/bar/bar'   | IllegalArgumentException  | "Invalid key in pod annotation -- Key can only contain exactly one '/' character"
    }

    @Unroll
    def 'should not sanitize k8s annotation value' () {
        given:
        def builder = new PodSpecBuilder()

        expect:
        builder.sanitize(ANNOTATION, PodSpecBuilder.MetaType.ANNOTATION) == EXPECTED

        where:
        ANNOTATION                      | EXPECTED
        ['foo':'value 1']               | ['foo':'value 1']
        ['foo':'foo.bar / *']           | ['foo':'foo.bar / *']
        ['foo':'value 2 \n value 3']    | ['foo':'value 2 \n value 3']
        ['foo':'value 3']               | ['foo':'value 3']
        ['foo':'x2345678901234567890123456789012345678901234567890123456789012345'] | ['foo':'x2345678901234567890123456789012345678901234567890123456789012345']
    }

    def 'should create job spec' () {

        when:
        def spec = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo', 'hello'])
                .buildAsJob()

        then:
        spec == [
            apiVersion: 'batch/v1',
            kind: 'Job',
            metadata: [name: 'foo', namespace: 'default'],
            spec: [
                backoffLimit: 0,
                template: [
                    metadata: [name: 'foo', namespace: 'default'],
                    spec: [
                        restartPolicy: 'Never',
                        containers: [[
                            name: 'foo',
                            image: 'busybox',
                            command: ['echo', 'hello'],
                        ]]
                    ]
                ]
            ]
        ]
    }

    def 'should create job spec with labels and annotations' () {

        when:
        def job = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo', 'hello'])
                .withLabel('app','someApp')
                .withLabel('runName','someName')
                .withLabel('version','3.8.1')
                .withAnnotation('anno1', 'val1')
                .withAnnotations([anno2: 'val2', anno3: 'val3'])
                .buildAsJob()

        def metadata = [
            name: 'foo',
            namespace: 'default',
            labels: [
                app: 'someApp',
                runName: 'someName',
                version: '3.8.1'
            ],
            annotations: [
                anno1: 'val1',
                anno2: 'val2',
                anno3: 'val3'
            ]
        ]

        then:
        job.metadata == metadata
        job.spec.template.metadata == metadata
    }

    def 'should create job spec with ttl seconds' () {
        when:
        def job = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo', 'hello'])
                .buildAsJob()
        then:
        !job.spec.ttlSecondsAfterFinished

        when:
        job = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo', 'hello'])
                .withPodOptions( new PodOptions(ttlSecondsAfterFinished: 60) )
                .buildAsJob()
        then:
        job.spec.ttlSecondsAfterFinished == 60
    }

}
