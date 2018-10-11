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
        options.getMountSecrets() == [] as Set
        options.getMountConfigMaps() == [] as Set
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

    def 'should create volume claims' () {
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
                [secret: 'secret/key', mountPath: '/etc/secret'],
                [config: 'data/key', mountPath: '/data/file.txt'],
                [volumeClaim: 'pvc', mountPath: '/mnt/claim'],
                [runAsUser: 500]
        ]

        def list2 = [
                [env: 'ALPHA', value: 'GAMMA'],
                [secret: 'foo/key', mountPath: '/a/aa'],
                [config: 'bar/key', mountPath: '/b/bb'],
                [volumeClaim: 'cvp', mountPath: '/c/cc'],

                [env: 'DELTA', value: 'LAMBDA'],
                [secret: 'x', mountPath: '/x'],
                [config: 'y', mountPath: '/y'],
                [volumeClaim: 'z', mountPath: '/z'],
        ]

        def list3 = [
                [env: 'HELLO', value: 'WORLD'],
                [secret: 'secret/key', mountPath: '/etc/secret'],
                [config: 'data/key', mountPath: '/data/file.txt'],
                [volumeClaim: 'pvc', mountPath: '/mnt/claim'],

                [env: 'DELTA', value: 'LAMBDA'],
                [secret: 'x', mountPath: '/x'],
                [config: 'y', mountPath: '/y'],
                [volumeClaim: 'z', mountPath: '/z'],
                [securityContext: [runAsUser: 1000, fsGroup: 200, allowPrivilegeEscalation: true]]

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

        opts.getMountSecrets() == [
            new PodMountSecret('secret/key', '/etc/secret'),
            new PodMountSecret('x', '/x')
        ] as Set

        opts.getMountConfigMaps() == [
                new PodMountConfig('data/key', '/data/file.txt'),
                new PodMountConfig('y', '/y'),
        ] as Set

        opts.getVolumeClaims() == [
                new PodVolumeClaim('pvc','/mnt/claim'),
                new PodVolumeClaim('z','/z'),
        ] as Set

        opts.securityContext.toSpec() == [runAsUser: 1000, fsGroup: 200, allowPrivilegeEscalation: true]
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
}
