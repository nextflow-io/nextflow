/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import nextflow.cli.CliOptions
import nextflow.cli.CmdKubeRun
import nextflow.cli.Launcher
import nextflow.k8s.client.ClientConfig
import nextflow.k8s.client.K8sClient
import nextflow.k8s.model.PodMountConfig
import nextflow.k8s.model.PodOptions
import nextflow.k8s.model.PodSpecBuilder
import nextflow.k8s.model.PodVolumeClaim
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sDriverLauncherTest extends Specification {

    def setup() {
        PodSpecBuilder.VOLUMES.set(0)
    }

    def 'should execute run' () {
        given:
        def NAME = 'nxf-foo'
        def NF_CONFIG = [process:[executor:'k8s']]
        def K8S_CONFIG = Mock(K8sConfig)
        def K8S_CLIENT = Mock(K8sClient)
        def driver = Spy(K8sDriverLauncher)

        when:
        driver.run(NAME, ['a','b','c'])
        then:
        1 * driver.makeConfig(NAME) >> NF_CONFIG
        1 * driver.makeK8sConfig(NF_CONFIG) >> K8S_CONFIG
        1 * driver.makeK8sClient(K8S_CONFIG) >> K8S_CLIENT
        1 * K8S_CONFIG.checkStorageAndPaths(K8S_CLIENT)
        1 * driver.createK8sConfigMap() >> null
        1 * driver.createK8sLauncherPod() >> null
        1 * driver.waitPodStart() >> null
        1 * driver.printK8sPodOutput() >> null

        driver.pipelineName == NAME
        driver.interactive == false
        driver.config == NF_CONFIG
        driver.k8sConfig == K8S_CONFIG
        driver.k8sClient == K8S_CLIENT
    }

    def 'should make k8s config' () {

        given:
        K8sConfig k8sConfig
        K8sDriverLauncher driver = Spy(K8sDriverLauncher)

        when:
        k8sConfig = driver.makeK8sConfig([:])
        then:
        k8sConfig == new K8sConfig()

        when:
        k8sConfig = driver.makeK8sConfig(k8s: [storageClaimName: 'foo', storageMountPath: '/mnt'])
        then:
        k8sConfig.getStorageClaimName() == 'foo'
        k8sConfig.getStorageMountPath() == '/mnt'

    }
 

    @Unroll
    def 'should get cmd cli' () {

        given:
        def l = new K8sDriverLauncher(cmd: cmd, pipelineName: 'foo')

        when:
        cmd.launcher = new Launcher(options: new CliOptions())
        then:
        l.getLaunchCli() == expected

        where:
        cmd                                         | expected
        new CmdKubeRun()                                | 'nextflow run foo'
        new CmdKubeRun(cacheable: false)                | 'nextflow run foo -cache false'
        new CmdKubeRun(resume: true)                    | 'nextflow run foo -resume true'
        new CmdKubeRun(poolSize: 10)                    | 'nextflow run foo -ps 10'
        new CmdKubeRun(pollInterval: 5)                 | 'nextflow run foo -pi 5'
        new CmdKubeRun(queueSize: 9)                    | 'nextflow run foo -qs 9'
        new CmdKubeRun(revision: 'xyz')                 | 'nextflow run foo -r xyz'
        new CmdKubeRun(latest: true)                    | 'nextflow run foo -latest true'
        new CmdKubeRun(withTrace: true)                 | 'nextflow run foo -with-trace true'
        new CmdKubeRun(withTimeline: true)              | 'nextflow run foo -with-timeline true'
        new CmdKubeRun(withDag: true)                   | 'nextflow run foo -with-dag true'
        new CmdKubeRun(dumpHashes: true)                | 'nextflow run foo -dump-hashes true'
        new CmdKubeRun(dumpChannels: 'lala')            | 'nextflow run foo -dump-channels lala'
        new CmdKubeRun(env: [XX:'hello', YY: 'world'])  | 'nextflow run foo -e.XX hello -e.YY world'
        new CmdKubeRun(process: [mem: '100',cpus:'2'])  | 'nextflow run foo -process.mem 100 -process.cpus 2'
        new CmdKubeRun(params: [alpha:'x', beta:'y'])   | 'nextflow run foo --alpha x --beta y'
        new CmdKubeRun(params: [alpha: '/path/*.txt'])  | 'nextflow run foo --alpha /path/\\*.txt'
    }

    def 'should set the run name' () {
        given:
        def cmd = new CmdKubeRun()
        cmd.launcher = new Launcher(options: new CliOptions())

        when:
        def l = new K8sDriverLauncher(cmd: cmd, pipelineName: 'foo', runName: 'bar')
        then:
        l.getLaunchCli() == 'nextflow run foo -name bar'
    }



    def 'should create launcher spec' () {

        given:
        def pod = Mock(PodOptions)
        pod.getVolumeClaims() >> [ new PodVolumeClaim('pvc-1', '/mnt/path/data') ]
        pod.getMountConfigMaps() >> [ new PodMountConfig('cfg-2', '/mnt/path/cfg') ]

        def k8s = Mock(K8sConfig)
        k8s.getNextflowImageName() >> 'the-image'
        k8s.getLaunchDir() >> '/the/user/dir'
        k8s.getWorkDir() >> '/the/work/dir'
        k8s.getProjectDir() >> '/the/project/dir'
        k8s.getPodOptions() >> pod

        def driver = Spy(K8sDriverLauncher)
        driver.runName = 'foo-boo'
        driver.k8sClient = new K8sClient(new ClientConfig(namespace: 'foo', serviceAccount: 'bar'))
        driver.k8sConfig = k8s

        when:
        def spec = driver.makeLauncherSpec()
        then:
        driver.getLaunchCli() >> 'nextflow run foo'

        spec == [apiVersion: 'v1',
                 kind: 'Pod',
                 metadata: [name:'foo-boo', namespace:'foo', labels:[app:'nextflow', runName:'foo-boo']],
                 spec: [restartPolicy:'Never',
                        containers:[
                                [name:'foo-boo',
                                 image:'the-image',
                                 command:['/bin/bash', '-c', "source /etc/nextflow/init.sh; nextflow run foo"],
                                 env:[
                                         [name:'NXF_WORK', value:'/the/work/dir'],
                                         [name:'NXF_ASSETS', value:'/the/project/dir'],
                                         [name:'NXF_EXECUTOR', value:'k8s']],
                                 volumeMounts:[
                                         [name:'vol-1', mountPath:'/mnt/path/data'],
                                         [name:'vol-2', mountPath:'/mnt/path/cfg']]]
                                ],
                        serviceAccountName:'bar',
                        volumes:[[name:'vol-1', persistentVolumeClaim:[claimName:'pvc-1']],
                                 [name:'vol-2', configMap:[name:'cfg-2'] ]]
                 ]
        ]

    }

    def 'should use user provided pod image' () {

        given:
        def pod = Mock(PodOptions)
        pod.getVolumeClaims() >> [ new PodVolumeClaim('pvc-1', '/mnt/path/data') ]
        pod.getMountConfigMaps() >> [ new PodMountConfig('cfg-2', '/mnt/path/cfg') ]

        def k8s = Mock(K8sConfig)
        k8s.getLaunchDir() >> '/the/user/dir'
        k8s.getWorkDir() >> '/the/work/dir'
        k8s.getProjectDir() >> '/the/project/dir'
        k8s.getPodOptions() >> pod

        def driver = Spy(K8sDriverLauncher)
        driver.runName = 'foo-boo'
        driver.k8sClient = new K8sClient(new ClientConfig(namespace: 'foo', serviceAccount: 'bar'))
        driver.k8sConfig = k8s
        driver.podImage = 'foo/bar'

        when:
        def spec = driver.makeLauncherSpec()
        then:
        driver.getLaunchCli() >> 'nextflow run foo'

        spec == [apiVersion: 'v1',
                 kind: 'Pod',
                 metadata: [name:'foo-boo', namespace:'foo', labels:[app:'nextflow', runName:'foo-boo']],
                 spec: [restartPolicy:'Never',
                        containers:[
                                [name:'foo-boo',
                                 image:'foo/bar',
                                 command:['/bin/bash', '-c', "source /etc/nextflow/init.sh; nextflow run foo"],
                                 env:[
                                         [name:'NXF_WORK', value:'/the/work/dir'],
                                         [name:'NXF_ASSETS', value:'/the/project/dir'],
                                         [name:'NXF_EXECUTOR', value:'k8s']],
                                 volumeMounts:[
                                         [name:'vol-1', mountPath:'/mnt/path/data'],
                                         [name:'vol-2', mountPath:'/mnt/path/cfg']]]
                        ],
                        serviceAccountName:'bar',
                        volumes:[[name:'vol-1', persistentVolumeClaim:[claimName:'pvc-1']],
                                 [name:'vol-2', configMap:[name:'cfg-2'] ]]
                 ]
        ]

    }

    def 'should create config map' () {

        given:
        def folder = Files.createTempDirectory('foo')

        def params = folder.resolve('params.json')
        params.text = 'bla-bla'
        def driver = Spy(K8sDriverLauncher)
        def NXF_CONFIG = [foo: 'bar']

        def SCM_FILE = folder.resolve('scm')
        SCM_FILE.text = "hello = 'world'\n"


        def EXPECTED = [:]
        EXPECTED['init.sh'] == ''

        def POD_OPTIONS = new PodOptions()

        def K8S_CONFIG = Mock(K8sConfig)
        K8S_CONFIG.getLaunchDir() >> '/launch/dir'
        K8S_CONFIG.getPodOptions() >> POD_OPTIONS

        when:
        driver.config = NXF_CONFIG
        driver.k8sConfig = K8S_CONFIG
        driver.cmd = new CmdKubeRun(paramsFile: params.toString())

        driver.createK8sConfigMap()
        then:
        1 * driver.getScmFile() >> SCM_FILE
        1 * driver.makeConfigMapName(_ as Map) >> 'nf-config-123'
        1 * driver.tryCreateConfigMap('nf-config-123', _ as Map) >> {  name, cfg ->
            assert cfg.'init.sh' == "mkdir -p '/launch/dir'; if [ -d '/launch/dir' ]; then cd '/launch/dir'; else echo 'Cannot create directory: /launch/dir'; exit 1; fi; [ -f /etc/nextflow/scm ] && ln -s /etc/nextflow/scm \$NXF_HOME/scm; [ -f /etc/nextflow/nextflow.config ] && cp /etc/nextflow/nextflow.config \$PWD/nextflow.config; "
            assert cfg.'nextflow.config' == "foo = 'bar'\n"
            assert cfg.'scm' == "hello = 'world'\n"
            assert cfg.'params.json' == 'bla-bla'
            return null
        }

        POD_OPTIONS.getMountConfigMaps() == [ new PodMountConfig('nf-config-123', '/etc/nextflow') ] as Set

        cleanup:
        folder?.deleteDir()
    }

    def 'should make config' () {
        given:
        Map config
        def driver = Spy(K8sDriverLauncher)
        def NAME = 'somePipelineName'
        def CFG_EMPTY = new ConfigObject()
        def CFG_WITH_MOUNTS = new ConfigObject()
        CFG_WITH_MOUNTS.k8s.storageClaimName = 'pvc'
        CFG_WITH_MOUNTS.k8s.storageMountPath = '/foo'

        when:
        driver.cmd = new CmdKubeRun()
        config = driver.makeConfig(NAME)
        then:
        1 *  driver.loadConfig(NAME) >> CFG_EMPTY
        config.process.executor == 'k8s'
        config.k8s.pod == null
        config.k8s.storageMountPath == null
        config.k8s.storageClaimName == null

        when:
        driver.cmd = new CmdKubeRun()
        config = driver.makeConfig(NAME)
        then:
        1 *  driver.loadConfig(NAME) >> CFG_WITH_MOUNTS
        config.process.executor == 'k8s'
        config.k8s.storageClaimName == 'pvc'
        config.k8s.storageMountPath == '/foo'
        and:
        new K8sConfig(config.k8s).getStorageClaimName() == 'pvc'
        new K8sConfig(config.k8s).getStorageMountPath() == '/foo'
        new K8sConfig(config.k8s).getPodOptions() == new PodOptions([ [volumeClaim:'pvc', mountPath: '/foo'] ])

        when:
        driver.cmd = new CmdKubeRun(volMounts: ['pvc-1:/this','pvc-2:/that'] )
        config = driver.makeConfig(NAME)
        then:
        1 *  driver.loadConfig(NAME) >> CFG_EMPTY
        config.process.executor == 'k8s'
        config.k8s.storageClaimName == 'pvc-1'
        config.k8s.storageMountPath == '/this'
        config.k8s.pod == [ [volumeClaim: 'pvc-2', mountPath: '/that'] ]
        and:
        new K8sConfig(config.k8s).getStorageClaimName() == 'pvc-1'
        new K8sConfig(config.k8s).getStorageMountPath() == '/this'
        new K8sConfig(config.k8s).getPodOptions() == new PodOptions([
                [volumeClaim:'pvc-1', mountPath: '/this'],
                [volumeClaim:'pvc-2', mountPath: '/that']
        ])


        when:
        driver.cmd = new CmdKubeRun(volMounts: ['xyz:/this'] )
        config = driver.makeConfig(NAME)
        then:
        1 *  driver.loadConfig(NAME) >> CFG_WITH_MOUNTS
        config.process.executor == 'k8s'
        config.k8s.storageClaimName == 'xyz'
        config.k8s.storageMountPath == '/this'
        config.k8s.pod == null
        and:
        and:
        new K8sConfig(config.k8s).getStorageClaimName() == 'xyz'
        new K8sConfig(config.k8s).getStorageMountPath() == '/this'
        new K8sConfig(config.k8s).getPodOptions() == new PodOptions([
                [volumeClaim:'xyz', mountPath: '/this']
        ])


        when:
        driver.cmd = new CmdKubeRun(volMounts: ['xyz', 'bar:/mnt/bar'] )
        config = driver.makeConfig(NAME)
        then:
        1 *  driver.loadConfig(NAME) >> CFG_WITH_MOUNTS
        config.process.executor == 'k8s'
        config.k8s.storageClaimName == 'xyz'
        config.k8s.storageMountPath == null
        config.k8s.pod == [ [volumeClaim: 'bar', mountPath: '/mnt/bar'] ]
        and:
        new K8sConfig(config.k8s).getStorageClaimName() == 'xyz'
        new K8sConfig(config.k8s).getStorageMountPath() == '/workspace'
        new K8sConfig(config.k8s).getPodOptions() == new PodOptions([
                [volumeClaim:'xyz', mountPath: '/workspace'],
                [volumeClaim:'bar', mountPath: '/mnt/bar']
        ])

    }

    def 'should make config - deprecated' () {

        given:
        Map config
        def driver = Spy(K8sDriverLauncher)
        def NAME = 'somePipelineName'
        def CFG_EMPTY = new ConfigObject()
        def CFG_WITH_MOUNTS = new ConfigObject()
        CFG_WITH_MOUNTS.k8s.volumeClaims = [ pvc: [mountPath:'/foo'] ]

        when:
        driver.cmd = new CmdKubeRun()
        config = driver.makeConfig(NAME)
        then:
        1 *  driver.loadConfig(NAME) >> CFG_EMPTY
        config.process.executor == 'k8s'

        when:
        driver.cmd = new CmdKubeRun()
        config = driver.makeConfig(NAME)
        then:
        1 *  driver.loadConfig(NAME) >> CFG_WITH_MOUNTS
        config.process.executor == 'k8s'
        config.k8s.storageClaimName == 'pvc'
        config.k8s.storageMountPath == '/foo'

        when:
        driver.cmd = new CmdKubeRun(volMounts: ['pvc-1:/this','pvc-2:/that'] )
        config = driver.makeConfig(NAME)
        then:
        1 *  driver.loadConfig(NAME) >> CFG_EMPTY
        config.process.executor == 'k8s'
        config.k8s.storageClaimName == 'pvc-1'
        config.k8s.storageMountPath == '/this'
        config.k8s.pod == [ [volumeClaim: 'pvc-2', mountPath: '/that'] ]


        when:
        driver.cmd = new CmdKubeRun(volMounts: ['xyz:/this'] )
        config = driver.makeConfig(NAME)
        then:
        1 *  driver.loadConfig(NAME) >> CFG_WITH_MOUNTS
        config.process.executor == 'k8s'
        config.k8s.storageClaimName == 'xyz'
        config.k8s.storageMountPath == '/this'
        config.k8s.pod == null
        and:
        new K8sConfig(config.k8s).getStorageClaimName() == 'xyz'
        new K8sConfig(config.k8s).getStorageMountPath() == '/this'
        new K8sConfig(config.k8s).getPodOptions() == new PodOptions([
                [volumeClaim:'xyz', mountPath: '/this']
        ])

    }
}
