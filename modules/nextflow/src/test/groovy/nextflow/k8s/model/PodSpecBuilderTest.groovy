/*
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

package nextflow.k8s.model

import nextflow.executor.res.AcceleratorResource
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


    def 'should set namespace, labels and annotations' () {

        when:
        def spec = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withWorkDir('/some/work/dir')
                .withCommand(['sh', '-c', 'echo hello'])
                .withNamespace('xyz')
                .withLabel('app','myApp')
                .withLabel('runName','something')
                .withAnnotation("anno1", "value1")
                .withAnnotations([anno2: "value2", anno3: "value3"])
                .build()

        then:
        spec ==  [ apiVersion: 'v1',
                   kind: 'Pod',
                   metadata: [
                           name:'foo',
                           namespace:'xyz',
                           labels: [
                                   app: 'myApp',
                                   runName: 'something'
                           ],
                           annotations: [
                                   anno1: "value1",
                                   anno2: "value2",
                                   anno3: "value3"
                           ]
                   ],
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
    
    def 'should truncate labels longer than 63 chars' () {

        when:
        def spec = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withWorkDir('/some/work/dir')
                .withCommand(['sh', '-c', 'echo hello'])
                .withNamespace('xyz')
                .withLabel('app','myApp')
                .withLabel('runName','something')
                .withLabel('tag','somethingreallylonggggggggggggggggggggggggggggggggggggggggggendEXTRABIT')
                .withLabels([tag2: 'somethingreallylonggggggggggggggggggggggggggggggggggggggggggendEXTRABIT', tag3: 'somethingreallylonggggggggggggggggggggggggggggggggggggggggggendEXTRABIT'])
                .withAnnotation("anno1", "value1")
                .withAnnotations([anno2: "value2", anno3: "value3"])
                .build()

        then:
        spec ==  [ apiVersion: 'v1',
                   kind: 'Pod',
                   metadata: [
                           name:'foo',
                           namespace:'xyz',
                           labels: [
                                   app: 'myApp',
                                   runName: 'something',
                                   tag: 'somethingreallylonggggggggggggggggggggggggggggggggggggggggggend',
                                   tag2: 'somethingreallylonggggggggggggggggggggggggggggggggggggggggggend',
                                   tag3: 'somethingreallylonggggggggggggggggggggggggggggggggggggggggggend'
                           ],
                           annotations: [
                                   anno1: "value1",
                                   anno2: "value2",
                                   anno3: "value3"
                           ]
                   ],
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
        def spec = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand('echo hello')
                .withWorkDir('/some/work/dir')
                .withEnv(PodEnv.value('ALPHA','hello'))
                .withEnv(PodEnv.value('DELTA', 'world'))
                .withCpus(8)
                .withAccelerator( new AcceleratorResource(request: 5, limit:10, type: 'foo.org') )
                .withMemory('100Gi')
                .build()

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
                                    resources:[
                                            requests: ['foo.org/gpu':5],
                                            limits:['foo.org/gpu':10, cpu:8, memory:'100Gi'] ]
                                   ]
                           ]
                   ]
        ]
    }

    def 'should get storage spec for volume claims' () {

        when:
        def spec = new PodSpecBuilder()
                    .withPodName('foo')
                    .withImageName('busybox')
                    .withWorkDir('/path')
                    .withCommand(['echo'])
                    .withVolumeClaim(new PodVolumeClaim('first','/work'))
                    .withVolumeClaim(new PodVolumeClaim('second', '/data', '/foo'))
                    .build()
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
                                            [name:'vol-2', mountPath:'/data', subPath: '/foo']] ]
                           ],
                           volumes:[
                                   [name:'vol-1', persistentVolumeClaim:[claimName:'first']],
                                   [name:'vol-2', persistentVolumeClaim:[claimName:'second']] ]
                   ]

        ]

    }

    def 'should only define one volume per persistentVolumeClaim' () {

        when:
        def spec = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withWorkDir('/path')
                .withCommand(['echo'])
                .withVolumeClaim(new PodVolumeClaim('first','/work'))
                .withVolumeClaim(new PodVolumeClaim('first','/work2', '/bar'))
                .withVolumeClaim(new PodVolumeClaim('second', '/data', '/foo'))
                .withVolumeClaim(new PodVolumeClaim('second', '/data2', '/fooz'))
                .build()
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
                                            [name:'vol-1', mountPath:'/work2', subPath: '/bar'],
                                            [name:'vol-2', mountPath:'/data', subPath: '/foo'],
                                            [name:'vol-2', mountPath:'/data2', subPath: '/fooz']]]
                           ],
                           volumes:[
                                   [name:'vol-1', persistentVolumeClaim:[claimName:'first']],
                                   [name:'vol-2', persistentVolumeClaim:[claimName:'second']] ]
                   ]

        ]

    }

    def 'should get config map mounts' () {

        when:
        def spec = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withWorkDir('/path')
                .withCommand(['echo'])
                .withConfigMap(new PodMountConfig(config: 'cfg1', mountPath: '/etc/config'))
                .withConfigMap(new PodMountConfig(config: 'data2', mountPath: '/data/path'))
                .build()
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

    def 'should consume env secrets' () {

        when:
        def spec = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo'])
                .withEnv( PodEnv.value('FOO','abc'))
                .withEnv( PodEnv.secret('VAR_X', 'delta/bar'))
                .withEnv( PodEnv.secret('VAR_Y', 'gamma'))
                .build()

        then:
        spec == [ apiVersion: 'v1',
                  kind: 'Pod',
                  metadata: [name:'foo', namespace:'default'],
                  spec: [
                          restartPolicy:'Never',
                          containers:[
                                  [name:'foo',
                                   image:'busybox',
                                   command: ['echo'],
                                   env:[   [name: 'FOO', value: 'abc'],
                                           [name:'VAR_X', valueFrom: [secretKeyRef: [name:'delta', key:'bar']]],
                                           [name:'VAR_Y', valueFrom: [secretKeyRef: [name:'gamma', key:'VAR_Y']]] ]
                                  ]
                          ]
                  ]

        ]
    }

    def 'should consume env configMap' () {

        when:
        def spec = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo'])
                .withEnv( PodEnv.value('FOO','abc'))
                .withEnv( PodEnv.config('VAR_X', 'data'))
                .withEnv( PodEnv.config('VAR_Y', 'omega/bar-2'))
                .build()

        then:
        spec == [ apiVersion: 'v1',
                  kind: 'Pod',
                  metadata: [name:'foo', namespace:'default'],
                  spec: [
                          restartPolicy:'Never',
                          containers:[
                                  [name:'foo',
                                   image:'busybox',
                                   command: ['echo'],
                                   env:[   [name: 'FOO', value: 'abc'],
                                           [name:'VAR_X', valueFrom: [configMapKeyRef: [name:'data', key:'VAR_X']]],
                                           [name:'VAR_Y', valueFrom: [configMapKeyRef: [name:'omega', key:'bar-2']]] ]
                                  ]
                          ]
                  ]

        ]
    }

    def 'should consume file secrets' () {

        when:
        def spec = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withCommand(['echo'])
                .withSecret(new PodMountSecret(secret: 'alpha', mountPath: '/this/and/that'))
                .withSecret(new PodMountSecret(secret: 'delta/foo', mountPath: '/etc/mnt/bar.txt'))
                .build()

        then:
        spec == [ apiVersion: 'v1',
                  kind: 'Pod',
                  metadata: [name:'foo', namespace:'default'],
                  spec: [
                          restartPolicy:'Never',
                          containers:[
                                  [name:'foo',
                                   image:'busybox',
                                   command: ['echo'],
                                   volumeMounts:[
                                           [name:'vol-1', mountPath:'/this/and/that'],
                                           [name:'vol-2', mountPath:'/etc/mnt']
                                   ]
                                  ]
                          ],
                          volumes:[
                                  [name:'vol-1', secret:[secretName: 'alpha']],
                                  [name:'vol-2', secret:[
                                          secretName: 'delta',
                                          items: [
                                                  [ key: 'foo', path:'bar.txt' ]
                                          ]
                                  ]]
                          ]
                  ]

        ]
    }

    def 'should get host path mounts' () {

        when:
        def spec = new PodSpecBuilder()
                .withPodName('foo')
                .withImageName('busybox')
                .withWorkDir('/path')
                .withCommand(['echo'])
                .withHostMount('/tmp','/scratch')
                .withHostMount('/host/data','/mnt/container')
                .build()

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


    def 'should create pod spec given pod options' () {

        given:
        def opts = Mock(PodOptions)
        def builder = new PodSpecBuilder(podName: 'foo', imageName: 'image', command: ['echo'], labels: [runName: 'crazy_john'], annotations: [evict: 'false'])

        when:
        def spec = builder.withPodOptions(opts).build()
        then:
        2 * opts.getImagePullPolicy() >> 'always'
        2 * opts.getImagePullSecret() >> 'myPullSecret'
        2 * opts.getVolumeClaims() >> [ new PodVolumeClaim('pvc1', '/work') ]
        2 * opts.getMountConfigMaps() >> [ new PodMountConfig('data', '/home/user') ]
        2 * opts.getMountSecrets() >> [ new PodMountSecret('blah', '/etc/secret.txt') ]
        2 * opts.getEnvVars() >> [ PodEnv.value('HELLO','WORLD') ]
        _ * opts.getLabels() >> [ALPHA: 'xxx', GAMMA: 'yyy']
        _ * opts.getAnnotations() >> [OMEGA:'zzz', SIGMA:'www']
        _ * opts.getSecurityContext() >> new PodSecurityContext(1000)
        _ * opts.getNodeSelector() >> new PodNodeSelector(gpu:true, queue: 'fast')

        spec == [
                apiVersion: 'v1',
                kind: 'Pod',
                metadata: [
                        name:'foo',
                        namespace:'default',
                        labels:[runName:'crazy_john', ALPHA:'xxx', GAMMA:'yyy'],
                        annotations: [evict: 'false', OMEGA:'zzz', SIGMA:'www']
                ],
                spec: [
                        restartPolicy:'Never',
                        securityContext: [ runAsUser: 1000 ],
                        imagePullSecrets: [[ name: 'myPullSecret' ]],
                        nodeSelector: [gpu: 'true', queue: 'fast'],
                        
                        containers:[
                                [name:'foo',
                                 image:'image',
                                        imagePullPolicy: 'always',
                                 command:['echo'],
                                 env:[[name:'HELLO', value:'WORLD']],
                                 volumeMounts:[
                                         [name:'vol-1', mountPath:'/work'],
                                         [name:'vol-2', mountPath:'/home/user'],
                                         [name:'vol-3', mountPath:'/etc/secret.txt']
                                 ],
                                ]
                        ],
                        volumes:[
                                [name:'vol-1', persistentVolumeClaim:[claimName:'pvc1']],
                                [name:'vol-2', configMap:[name:'data']],
                                [name:'vol-3', secret:[secretName:'blah']] ]
                ]

        ]

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
        res = builder.addAcceleratorResources(new AcceleratorResource(request: 5, type:'foo.org'), [limits: [cpus: 2]])
        then:
        res.requests == ['foo.org/gpu': 5]
        res.limits == [cpus:2]

        when:
        res = builder.addAcceleratorResources(new AcceleratorResource(request: 5, limit: 10, type:'foo.org'), [limits: [cpus: 2]])
        then:
        res.requests == ['foo.org/gpu': 5]
        res.limits == [cpus:2, 'foo.org/gpu': 10]

    }


    @Unroll
    def 'should sanitize k8s label: #label' () {
        given:
        def builder = new PodSpecBuilder()

        expect:
        builder.sanitize0('foo',label, 'label') == str

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
        'x2345678901234567890123456789012345678901234567890123456789012345' | 'x23456789012345678901234567890123456789012345678901234567890123'
    }
}
