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

package nextflow.k8s

import nextflow.BuildInfo
import nextflow.k8s.client.ClientConfig
import nextflow.k8s.model.PodEnv
import nextflow.k8s.model.PodSecurityContext
import nextflow.k8s.model.PodVolumeClaim
import nextflow.util.Duration
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sConfigTest extends Specification {

    def 'should create config object' () {

        when:
        def cfg = new K8sConfig()
        then:
        cfg.getNamespace() == null
        cfg.getServiceAccount() == null
        !cfg.getDebug().getYaml()

        when:
        cfg = new K8sConfig( namespace:'foo', serviceAccount: 'bar', debug: [yaml: true] )
        then:
        cfg.getNamespace() == 'foo'
        cfg.getServiceAccount() == 'bar'
        cfg.getDebug().getYaml()
        cfg.debug.yaml

    }

    def 'should set cleanup' () {
        given:
        K8sConfig cfg

        when:
        cfg = new K8sConfig()
        then: 'it should return true when missing value'
        cfg.getCleanup()

        when:
        cfg = new K8sConfig()
        then: 'it should return false specified as default'
        !cfg.getCleanup(false)

        when:
        cfg = new K8sConfig(cleanup:false)
        then: 'it should return false'
        !cfg.getCleanup()

        when:
        cfg = new K8sConfig(cleanup:'false')
        then: 'it should return false'
        !cfg.getCleanup()

        when:
        cfg = new K8sConfig(cleanup:true)
        then: 'it should return true'
        cfg.getCleanup()

        when:
        cfg = new K8sConfig(cleanup:'true')
        then: 'it should return true'
        cfg.getCleanup()

        when:
        cfg = new K8sConfig(cleanup:true)
        then: 'the default value should be ignored'
        cfg.getCleanup(false)

        when:
        cfg = new K8sConfig(cleanup:'true')
        then: 'the default value should be ignored'
        cfg.getCleanup(false)
    }

    def 'should create config with storage claims' () {

        when:
        def cfg = new K8sConfig(storageClaimName: 'pvc-1')
        then:
        cfg.getStorageClaimName() == 'pvc-1'
        cfg.getStorageMountPath() == '/workspace'
        cfg.getPodOptions().getVolumeClaims() == [ new PodVolumeClaim('pvc-1', '/workspace') ] as Set

        when:
        cfg = new K8sConfig([
                storageClaimName: 'pvc-2',
                storageMountPath: '/data',
                pod: [  [volumeClaim:'foo', mountPath: '/here'],
                        [volumeClaim: 'bar', mountPath: '/there']] ])
        then:
        cfg.getStorageClaimName() == 'pvc-2'
        cfg.getStorageMountPath() == '/data'
        cfg.getPodOptions().getVolumeClaims() == [
                new PodVolumeClaim('pvc-2', '/data'),
                new PodVolumeClaim('foo', '/here'),
                new PodVolumeClaim('bar', '/there')
        ] as Set


        when:
        cfg = new K8sConfig(storageClaimName: 'pvc-3', storageMountPath: '/some/path', storageSubPath: '/bar')
        then:
        cfg.getStorageClaimName() == 'pvc-3'
        cfg.getStorageMountPath() == '/some/path'
        cfg.getStorageSubPath() == '/bar'
        cfg.getPodOptions().getVolumeClaims() == [ new PodVolumeClaim('pvc-3', '/some/path', '/bar') ] as Set

    }

    def 'should set device plugin' () {
        when:
        def cfg = new K8sConfig([:])
        then:
        cfg.fuseDevicePlugin() == ['nextflow.io/fuse':1]

        when:
        cfg = new K8sConfig([fuseDevicePlugin:['foo/fuse':10]])
        then:
        cfg.fuseDevicePlugin() == ['foo/fuse':10]
    }

    def 'should create client config' () {

        given:
        def CONFIG = [namespace: 'this', serviceAccount: 'that', client: [server: 'http://foo']]

        when:
        def config = new K8sConfig(CONFIG)
        def client = config.getClient()
        then:
        client.server == 'http://foo'
        client.namespace == 'this'
        client.serviceAccount == 'that'
        client.httpConnectTimeout == null // testing default null
        client.httpReadTimeout == null // testing default null
        client.maxErrorRetry == 4

    }

    def 'should set maxErrorRetry' () {
        given:
        def CONFIG = [maxErrorRetry: 10, namespace: 'this', serviceAccount: 'that', client: [server: 'http://foo']]

        when:
        def config = new K8sConfig(CONFIG)
        def client = config.getClient()
        then:
        client.maxErrorRetry == 10
    }

    def 'should create client config with http request timeouts' () {

        given:
        def CONFIG = [
                namespace: 'this',
                serviceAccount: 'that',
                client: [server: 'http://foo'],
                httpReadTimeout: '20s',
                httpConnectTimeout: '25s' ]

        when:
        def config = new K8sConfig(CONFIG)
        def client = config.getClient()
        then:
        client.server == 'http://foo'
        client.namespace == 'this'
        client.serviceAccount == 'that'
        client.httpConnectTimeout == Duration.of('25s')
        client.httpReadTimeout == Duration.of('20s')

    }

    @Unroll
    def 'should create client config with discovery' () {

        given:
        def CONFIG = [context: CONTEXT, namespace: NAMESPACE, serviceAccount: SERVICE_ACCOUNT]
        K8sConfig config = Spy(K8sConfig, constructorArgs: [ CONFIG ])

        when:
        def client = config.getClient()
        then:
        1 * config.clientDiscovery(CONTEXT, NAMESPACE, SERVICE_ACCOUNT) >> new ClientConfig(namespace: NAMESPACE, server: SERVER)
        and:
        client.server == SERVER
        client.namespace == NAMESPACE ?: 'default'
        client.serviceAccount == SERVICE_ACCOUNT ?: 'default'

        where:
        CONTEXT     | SERVER    | NAMESPACE | SERVICE_ACCOUNT
        'foo'       | 'host.com'| null      | null
        'bar'       | 'this.com'| 'ns1'     | 'sa2'

    }

    def 'should get nextflow image name' () {

        when:
        def cfg = new K8sConfig()
        then:
        cfg.getNextflowImageName() ==  "nextflow/nextflow:${BuildInfo.version}"

        when:
        cfg = new K8sConfig(nextflow: [image: 'foo/bar:1.0'])
        then:
        cfg.getNextflowImageName() == 'foo/bar:1.0'

    }

    def 'should get autoMountHostPaths' () {

        when:
        def cfg = new K8sConfig()
        then:
        !cfg.getAutoMountHostPaths()

        when:
        cfg = new K8sConfig(autoMountHostPaths: true)
        then:
        cfg.getAutoMountHostPaths()

        when:
        cfg = new K8sConfig(autoMountHostPaths: 'true')
        then:
        cfg.getAutoMountHostPaths()

        when:
        cfg = new K8sConfig(autoMountHostPaths: false)
        then:
        !cfg.getAutoMountHostPaths()

        when:
        cfg = new K8sConfig(autoMountHostPaths: 'false')
        then:
        !cfg.getAutoMountHostPaths()
    }


    def 'should get podOptions' () {

        when:
        def cfg = new K8sConfig()
        def opts = cfg.getPodOptions()
        then:
        opts.envVars == [] as Set
        opts.mountSecrets == [] as Set
        opts.mountConfigMaps == [] as Set
        opts.volumeClaims == [] as Set


        when:
        opts = new K8sConfig(pod: [ [pullPolicy: 'Always'], [env: 'HELLO', value: 'WORLD'] ]).getPodOptions()
        then:
        opts.getImagePullPolicy() == 'Always'
        opts.getEnvVars() == [ PodEnv.value('HELLO','WORLD') ] as Set
    }

    def 'should return user name' () {

        when:
        def cfg = new K8sConfig()
        then:
        cfg.getUserName() == System.properties.get('user.name')

        when:
        cfg = new K8sConfig(userName: 'foo')
        then:
        cfg.getUserName() == 'foo'
    }

    def 'should return user dir' () {
        when:
        def cfg = new K8sConfig()
        then:
        cfg.getLaunchDir() == '/workspace/' + System.properties.get('user.name')

        when:
        cfg = new K8sConfig(storageMountPath: '/this/path', userName: 'foo')
        then:
        cfg.getLaunchDir() == '/this/path/foo'

        when:
        cfg = new K8sConfig(storageMountPath: '/this/path', userName: 'foo', launchDir: '/my/path')
        then:
        cfg.getLaunchDir() == '/my/path'

    }

    def 'should return work dir' () {
        when:
        def cfg = new K8sConfig()
        then:
        cfg.getWorkDir() == "/workspace/${System.properties.get('user.name')}/work"

        when:
        cfg = new K8sConfig(launchDir: '/my/dir')
        then:
        cfg.getWorkDir() == "/my/dir/work"

        when:
        cfg = new K8sConfig(launchDir: '/my/dir', workDir: '/the/wor/dir')
        then:
        cfg.getWorkDir() == "/the/wor/dir"
    }

    def 'should return project dir' () {
        when:
        def cfg = new K8sConfig()
        then:
        cfg.getProjectDir() == '/workspace/projects'

        when:
        cfg = new K8sConfig(storageMountPath: '/foo')
        then:
        cfg.getProjectDir() == '/foo/projects'

        when:
        cfg = new K8sConfig(storageMountPath: '/foo', projectDir: '/my/project/dir')
        then:
        cfg.getProjectDir() == '/my/project/dir'
    }

    def 'should return storage dir' () {

        when:
        def cfg = new K8sConfig()
        then:
        cfg.getStorageMountPath() == '/workspace'

        when:
        cfg = new K8sConfig(storageMountPath: '/mnt/there')
        then:
        cfg.getStorageMountPath() == '/mnt/there'

    }

    def 'should return compute resource type' () {

        when:
        def cfg = new K8sConfig()
        then:
        !cfg.useJobResource()

        when:
        cfg = new K8sConfig(computeResourceType: 'Job')
        then:
        cfg.useJobResource()

    }

    def 'should return storage claim name' () {
        when:
        def cfg = new K8sConfig()
        then:
        cfg.getStorageClaimName() == null

        when:
        cfg = new K8sConfig(storageClaimName: 'xxx')
        then:
        cfg.getStorageClaimName() == 'xxx'
    }

    def 'should create k8s config with one volume claim' () {

        when:
        def cfg = new K8sConfig( pod: [runAsUser: 1000] )
        then:
        cfg.getPodOptions().getSecurityContext() == new PodSecurityContext(1000)
        cfg.getPodOptions().getVolumeClaims().size() == 0

        when:
        cfg = new K8sConfig( pod: [volumeClaim: 'nf-0001', mountPath: '/workspace'] )
        then:
        cfg.getPodOptions().getSecurityContext() == null
        cfg.getPodOptions().getVolumeClaims() == [new PodVolumeClaim('nf-0001', '/workspace')] as Set


        when:
        cfg = new K8sConfig( pod: [
                [runAsUser: 1000],
                [volumeClaim: 'nf-0001', mountPath: '/workspace'],
                [volumeClaim: 'nf-0002', mountPath: '/data', subPath: '/home']
        ])
        then:
        cfg.getPodOptions().getSecurityContext() == new PodSecurityContext(1000)
        cfg.getPodOptions().getVolumeClaims() == [
                    new PodVolumeClaim('nf-0001', '/workspace'),
                    new PodVolumeClaim('nf-0002', '/data', '/home')
        ] as Set
        
    }


    def 'should set the sec context'( ) {

        given:
        def ctx = [runAsUser: 500, fsGroup: 200, allowPrivilegeEscalation: true, seLinuxOptions: [level: "s0:c123,c456"]]

        when:
        def cfg = new K8sConfig( runAsUser: 500 )
        then:
        cfg.getPodOptions().getSecurityContext() == new PodSecurityContext(500)

        when:
        cfg = new K8sConfig( securityContext: ctx )
        then:
        cfg.getPodOptions().getSecurityContext() == new PodSecurityContext(ctx)

    }

    def 'should set env and sec context' () {
        given:
        def ctx = [
                [env: 'FUSION_BUCKETS', value: 's3://nextflow-ci'],
                [securityContext: [privileged: true]]]

        when:
        def cfg = new K8sConfig( pod: ctx )
        then:
        cfg.getPodOptions().getEnvVars().first() == PodEnv.value('FUSION_BUCKETS', 's3://nextflow-ci')
        cfg.getPodOptions().getSecurityContext().toSpec() == [privileged:true]

    }

    def 'should set the image pull policy' () {
        when:
        def cfg = new K8sConfig( pullPolicy: 'always' )
        then:
        cfg.getPodOptions().getImagePullPolicy() == 'always'
    }

    def 'should set preserve entrypoint setting'( ) {

        when:
        def cfg = new K8sConfig([:])
        then:
        !cfg.entrypointOverride()

        when:
        cfg = new K8sConfig( entrypointOverride: true )
        then:
        cfg.entrypointOverride()

    }

    def 'should set debug.yaml' () {
        when:
        def cfg = new K8sConfig( debug: [yaml: 'true'] )
        then:
        cfg.getDebug().getYaml()

        when:
        cfg = new K8sConfig( debug: [yaml: true] )
        then:
        cfg.getDebug().getYaml()

        when:
        cfg = new K8sConfig( debug: [yaml: 'false'] )
        then:
        !cfg.getDebug().getYaml()

        when:
        cfg = new K8sConfig( debug: [yaml: false] )
        then:
        !cfg.getDebug().getYaml()

        when:
        cfg = new K8sConfig( debug: null )
        then:
        !cfg.getDebug().getYaml()

        when:
        cfg = new K8sConfig( debug: [:] )
        then:
        !cfg.getDebug().getYaml()

        when:
        cfg = new K8sConfig( null )
        then:
        !cfg.getDebug().getYaml()
    }
  
    def 'should set fetchNodeName' () {
        when:
        def cfg = new K8sConfig( fetchNodeName: 'true' )
        then:
        cfg.fetchNodeName() == true

        when:
        cfg = new K8sConfig( fetchNodeName: true )
        then:
        cfg.fetchNodeName() == true

        when:
        cfg = new K8sConfig( fetchNodeName: 'false' )
        then:
        cfg.fetchNodeName() == false

        when:
        cfg = new K8sConfig( fetchNodeName: false )
        then:
        cfg.fetchNodeName() == false

        when:
        cfg = new K8sConfig()
        then:
        cfg.fetchNodeName() == false
    }
}
