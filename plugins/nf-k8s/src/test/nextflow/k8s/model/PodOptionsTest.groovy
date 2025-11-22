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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PodOptionsTest extends Specification {


    def 'should create empty options' () {

        when:
        def options = new PodOptions(null)
        then:
        options.getEnvVars() == [] as Set
        options.getMountConfigMaps() == [] as Set
        options.getMountCsiEphemerals() == [] as Set
        options.getMountEmptyDirs() == [] as Set
        options.getMountSecrets() == [] as Set
        options.getAutomountServiceAccountToken() == true
    }

    def 'should set pullPolicy' () {

        when:
        def options = new PodOptions()
        then:
        options.getImagePullPolicy() == null
        
        when:
        options = new PodOptions([ [pullPolicy:'Always'] ])
        then:
        options.getImagePullPolicy() == 'Always'

        when:
        options = new PodOptions([ [imagePullPolicy:'latest'] ])
        then:
        options.getImagePullPolicy() == 'latest'
    }

    def 'should set imagePullSecret' () {

        when:
        def options = new PodOptions()
        then:
        options.imagePullSecret == null

        when:
        options = new PodOptions([ [imagePullSecret:'foo'] ])
        then:
        options.imagePullSecret == 'foo'

        when:
        options = new PodOptions([ [imagePullSecrets:'bar'] ])
        then:
        options.imagePullSecret == 'bar'
    }

    def 'should return config mounts' () {

        given:
        def options = [
                [mountPath: '/this/path1.txt', config: 'name/key1'],
                [mountPath: '/this/path2.txt', config: 'name/key2'],
                [mountPath: '/this/path2.txt', config: 'name/key2'],  // <-- identical entry are ignored
                [mountPath: '/this/path2.txt', secret: 'name/secret'],
                [env: 'FOO', config: '/name/foo']
        ]

        when:
        def configs = new PodOptions(options).getMountConfigMaps()
        then:
        configs.size() == 2
        configs == [
                new PodMountConfig(mountPath: '/this/path1.txt', config: 'name/key1'),
                new PodMountConfig(mountPath: '/this/path2.txt', config: 'name/key2')
        ] as Set
    }

    def 'should return csi ephemeral mounts' () {

        given:
        def options = [
            [
                mountPath: '/data',
                csi: [
                    driver: 'inline.storage.kubernetes.io',
                    volumeAttributes: [foo: 'bar']
                ]
            ]
        ]

        when:
        def csiEphemerals = new PodOptions(options).getMountCsiEphemerals()
        then:
        csiEphemerals == [
            new PodMountCsiEphemeral(mountPath: '/data', csi: options[0].csi)
        ] as Set
    }

    def 'should return emptyDir mounts' () {

        given:
        def options = [
                [mountPath: '/scratch1', emptyDir: [medium: 'Memory']],
                [mountPath: '/scratch2', emptyDir: [medium: 'Disk']]
        ]

        when:
        def emptyDirs = new PodOptions(options).getMountEmptyDirs()
        then:
        emptyDirs.size() == 2
        emptyDirs == [
                new PodMountEmptyDir(options[0]),
                new PodMountEmptyDir(options[1]) ] as Set
    }

    def 'should return secret mounts' () {

        given:
        def options = [
                [mountPath: '/this/path1.txt', config: 'name/key1'],
                [mountPath: '/this/alpha.txt', secret: 'name/secret1'],
                [mountPath: '/this/beta.txt', secret: 'name/secret2'],
                [mountPath: '/this/beta.txt', secret: 'name/secret2'],
                [env: 'FOO', config: '/name/foo']
        ]

        when:
        def secrets = new PodOptions(options).getMountSecrets()
        then:
        secrets.size() == 2
        secrets == [
                new PodMountSecret(mountPath: '/this/alpha.txt', secret: 'name/secret1'),
                new PodMountSecret(mountPath: '/this/beta.txt', secret: 'name/secret2'),
        ] as Set
    }


    def 'should return env definitions' () {

        given:
        def options = [
                [mountPath: '/this/path1.txt', config: 'name/key1'],
                [mountPath: '/this/alpha.txt', secret: 'name/secret1'],
                [mountPath: '/this/beta.txt', secret: 'name/secret2'],
                [env: 'FOO', config: '/name/foo'],
                [env: 'FOO', config: '/name/foo'],
                [env: 'BAR', config: '/name/BAR'],
                [env: 'ALPHA', value: 'aaa'],
                [env: 'ALPHA', value: 'aaa'],
                [env: 'BETA',  value: 'bbb'],
                [env: 'PASSWORD', secret:'name/key'],
                [env: 'PASSWORD', secret:'name/key'],
        ]

        when:
        def env = new PodOptions(options).getEnvVars()
        then:
        env.size() == 5
        env == [
                PodEnv.config('FOO', '/name/foo'),
                PodEnv.config('BAR', '/name/BAR'),
                PodEnv.value('ALPHA', 'aaa'),
                PodEnv.value('BETA',  'bbb'),
                PodEnv.secret('PASSWORD', 'name/key'),
        ] as Set 

    }

    def 'should create persistent volume claims' () {
        given:
        def options = [
                [volumeClaim:'pvc1', mountPath: '/this/path'],
                [volumeClaim:'pvc2', mountPath: '/that/path'],
                [volumeClaim:'pvc3', mountPath: '/some/data', subPath: '/foo']
        ]

        when:
        def claims = new PodOptions(options).getVolumeClaims()
        then:
        claims.size() == 3
        claims == [
                new PodVolumeClaim('pvc1', '/this/path'),
                new PodVolumeClaim('pvc2', '/that/path'),
                new PodVolumeClaim('pvc3', '/some/data', '/foo')
        ] as Set

    }

    def 'should create host path' () {
        given:
        def options = [
            [hostPath: '/host/one', mountPath: '/pod/1'],
            [hostPath: '/host/two', mountPath: '/pod/2']
        ]
        when:
        def mounts = new PodOptions(options).getMountHostPaths()

        then:
        mounts == [
            new PodHostMount('/host/one', '/pod/1'),
            new PodHostMount('/host/two', '/pod/2')
        ] as Set

    }

    def 'should not create env' () {
        when:
        new PodOptions([ [env:'FOO'] ])
        then:
        thrown(IllegalArgumentException)

        when:
        new PodOptions([ [secret:'FOO'] ])
        then:
        thrown(IllegalArgumentException)

        when:
        new PodOptions([ [config:'FOO'] ])
        then:
        thrown(IllegalArgumentException)

        when:
        new PodOptions([ [volumeClaim:'FOO'] ])
        then:
        thrown(IllegalArgumentException)
    }

    def 'should merge podOptions' () {

        given:
        def list1 = [
                [env: 'HELLO', value: 'WORLD'],
                [config: 'data/key', mountPath: '/data/file.txt'],
                [secret: 'secret/key', mountPath: '/etc/secret'],
                [volumeClaim: 'pvc', mountPath: '/mnt/claim'],
                [runAsUser: 500]
        ]

        def list2 = [
                [env: 'ALPHA', value: 'GAMMA'],
                [config: 'bar/key', mountPath: '/b/bb'],
                [secret: 'foo/key', mountPath: '/a/aa'],
                [volumeClaim: 'cvp', mountPath: '/c/cc'],

                [env: 'DELTA', value: 'LAMBDA'],
                [config: 'y', mountPath: '/y'],
                [secret: 'x', mountPath: '/x'],
                [volumeClaim: 'z', mountPath: '/z'],
        ]

        def list3 = [
                [env: 'HELLO', value: 'WORLD'],
                [config: 'data/key', mountPath: '/data/file.txt'],
                [secret: 'secret/key', mountPath: '/etc/secret'],
                [volumeClaim: 'pvc', mountPath: '/mnt/claim'],

                [env: 'DELTA', value: 'LAMBDA'],
                [config: 'y', mountPath: '/y'],
                [secret: 'x', mountPath: '/x'],
                [volumeClaim: 'z', mountPath: '/z'],

                [csi: [driver: 'inline.storage.kubernetes.io'], mountPath: '/data'],
                [emptyDir: [:], mountPath: '/scratch1'],
                [securityContext: [runAsUser: 1000, fsGroup: 200, allowPrivilegeEscalation: true]],
                [nodeSelector: 'foo=X, bar=Y'],
                [automountServiceAccountToken: false],
                [priorityClassName: 'high-priority']
        ]

        PodOptions opts

        when:
        opts = new PodOptions() + new PodOptions()
        then:
        opts == new PodOptions()

        when:
        opts = new PodOptions(list1) + new PodOptions()
        then:
        opts == new PodOptions(list1)
        opts.securityContext.toSpec() == [runAsUser:500]

        when:
        opts = new PodOptions() + new PodOptions(list1)
        then:
        opts == new PodOptions(list1)
        opts.securityContext.toSpec() == [runAsUser:500]

        when:
        opts = new PodOptions(list1) + new PodOptions(list1)
        then:
        opts == new PodOptions(list1)
        opts.securityContext.toSpec() == [runAsUser:500]

        when:
        opts = new PodOptions(list1) + new PodOptions(list2)
        then:
        opts == new PodOptions(list1 + list2)
        opts.securityContext.toSpec() == [runAsUser:500]

        when:
        opts = new PodOptions(list1) + new PodOptions(list3)
        then:
        opts.getEnvVars() == [
            PodEnv.value('HELLO','WORLD'),
            PodEnv.value('DELTA','LAMBDA')
        ] as Set

        opts.getMountConfigMaps() == [
            new PodMountConfig('data/key', '/data/file.txt'),
            new PodMountConfig('y', '/y'),
        ] as Set

        opts.getMountCsiEphemerals() == [
            new PodMountCsiEphemeral([driver: 'inline.storage.kubernetes.io'], '/data')
        ] as Set

        opts.getMountEmptyDirs() == [
            new PodMountEmptyDir([:], '/scratch1'),
        ] as Set

        opts.getMountSecrets() == [
            new PodMountSecret('secret/key', '/etc/secret'),
            new PodMountSecret('x', '/x')
        ] as Set

        opts.getVolumeClaims() == [
            new PodVolumeClaim('pvc','/mnt/claim'),
            new PodVolumeClaim('z','/z'),
        ] as Set

        opts.securityContext.toSpec() == [runAsUser: 1000, fsGroup: 200, allowPrivilegeEscalation: true]

        opts.nodeSelector.toSpec() == [foo: 'X', bar: "Y"]

        opts.getAutomountServiceAccountToken() == false

        opts.getPriorityClassName() == 'high-priority'
    }

    def 'should copy image pull policy' (){
        given:
        def data = [
            [imagePullPolicy : 'FOO']
        ]

        when:
        def opts = new PodOptions() + new PodOptions(data)
        then:
        opts.imagePullPolicy == 'FOO'

        when:
        opts = new PodOptions(data) + new PodOptions()
        then:
        opts.imagePullPolicy == 'FOO'
    }

    def 'should copy image pull secret' (){
        given:
        def data = [
                [imagePullSecret : 'BAR']
        ]

        when:
        def opts = new PodOptions() + new PodOptions(data)
        then:
        opts.imagePullSecret == 'BAR'

        when:
        opts = new PodOptions(data) + new PodOptions()
        then:
        opts.imagePullSecret == 'BAR'
    }

    def 'should copy pod labels' (){
        given:
        def data = [
                [label: "LABEL", value: 'VALUE']
        ]

        when:
        def opts = new PodOptions() + new PodOptions(data)
        then:
        opts.labels == ["LABEL": "VALUE"]

        when:
        opts = new PodOptions(data) + new PodOptions()
        then:
        opts.labels == ["LABEL": "VALUE"]

        when:
        opts = new PodOptions([[label:"FOO", value:'one']]) + new PodOptions([[label:"BAR", value:'two']])
        then:
        opts.labels == [FOO: 'one', BAR: 'two']
    }

    def 'should copy host paths' (){
        given:
        def data = [
            [hostPath: "/foo", mountPath: '/one']
        ]

        when:
        def opts = new PodOptions() + new PodOptions(data)
        then:
        opts.getMountHostPaths() == [new PodHostMount('/foo', '/one')] as Set

        when:
        opts = new PodOptions(data) + new PodOptions()
        then:
        opts.getMountHostPaths() == [new PodHostMount('/foo', '/one')] as Set

        when:
        opts = new PodOptions([[hostPath:"/foo", mountPath: '/one']]) + new PodOptions([[hostPath:"/bar", mountPath: '/two']])
        then:
        opts.getMountHostPaths() == [
            new PodHostMount('/foo','/one'),
            new PodHostMount('/bar','/two')
        ] as Set
    }

    def 'should create pod labels' () {

        given:
        def options = [
                [label: 'ALPHA', value: 'aaa'],
                [label: 'DELTA', value: 'bbb'],
                [label: 'DELTA', value: 'ddd']
        ]
        
        when:
        def opts = new PodOptions(options)
        then:
        opts.labels.size() == 2
        opts.labels == [ALPHA: 'aaa', DELTA: 'ddd']

    }

    def 'should copy pod annotations' (){
        given:
        def data = [
                [annotation: "ANNOTATION", value: 'VALUE']
        ]

        when:
        def opts = new PodOptions() + new PodOptions(data)
        then:
        opts.annotations == ["ANNOTATION": "VALUE"]

        when:
        opts = new PodOptions(data) + new PodOptions()
        then:
        opts.annotations == ["ANNOTATION": "VALUE"]

        when:
        opts = new PodOptions([[annotation:"FOO", value:'one']]) + new PodOptions([[annotation:"BAR", value:'two']])
        then:
        opts.annotations == [FOO: 'one', BAR: 'two']
    }

    def 'should create pod annotations' () {

        given:
        def options = [
                [annotation: 'ALPHA', value: 'aaa'],
                [annotation: 'DELTA', value: 'bbb'],
                [annotation: 'DELTA', value: 'ddd']
        ]

        when:
        def opts = new PodOptions(options)
        then:
        opts.annotations.size() == 2
        opts.annotations == [ALPHA: 'aaa', DELTA: 'ddd']

    }

    def 'should create user security context' () {
        when:
        def opts = new PodOptions([ [runAsUser: 1000] ])
        then:
        opts.getSecurityContext() == new PodSecurityContext(1000)

        when:
        opts = new PodOptions([ [runAsUser: 'foo'] ])
        then:
        opts.getSecurityContext() == new PodSecurityContext('foo')

        when:
        opts = new PodOptions([ [runAsUser: 'foo'] ])
        then:
        opts.getSecurityContext() != new PodSecurityContext('bar')

        when:
        def ctx = [runAsUser: 500, fsGroup: 200, allowPrivilegeEscalation: true, seLinuxOptions: [level: "s0:c123,c456"]]
        def expected = new PodSecurityContext(ctx)
        opts = new PodOptions([ [securityContext: ctx] ])
        then:
        opts.getSecurityContext() == expected
        opts.getSecurityContext().toSpec() == ctx
    }

    def 'should create pod node selector' () {
        when:
        def opts = new PodOptions([ [nodeSelector: 'foo=1, bar=true, baz=Z'] ])
        then:
        opts.nodeSelector.toSpec() == [foo: '1', bar: 'true', baz: 'Z']

    }

    def 'should set pod automount service token' () {
        when:
        def opts = new PodOptions([[automountServiceAccountToken: false]])
        then:
        opts.getAutomountServiceAccountToken() == false
    }

    def 'should set pod priority class name' () {
        when:
        def opts = new PodOptions([[priorityClassName: 'high-priority']])
        then:
        opts.getPriorityClassName() == 'high-priority'
    }

    def 'should set pod privileged' () {
        when:
        def opts = new PodOptions([:])
        then:
        !opts.getPrivileged()

        when:
        opts = new PodOptions([[privileged: true]])
        then:
        opts.getPrivileged()
    }

    def 'should set pod schedulerName' () {
        when:
        def opts = new PodOptions()
        then:
        opts.getSchedulerName() == null

        when:
        opts = new PodOptions([ [schedulerName:'my-scheduler'] ])
        then:
        opts.getSchedulerName() == 'my-scheduler'
    }
}
