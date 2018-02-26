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

import java.nio.file.Files

import nextflow.cli.CliOptions
import nextflow.cli.CmdKubeRun
import nextflow.cli.Launcher
import nextflow.k8s.client.ClientConfig
import nextflow.k8s.client.K8sClient
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sDriverLauncherTest extends Specification {

    def setup() {
        K8sHelper.VOLUMES.set(0)
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
        new CmdKubeRun(profile: 'ciao')                 | 'nextflow run foo -profile ciao'
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


    def 'should create config' () {

        given:
        Map cfg
        K8sClient client
        def driver = Spy(K8sDriverLauncher)
        def CLIENT_CFG = [server: 'foo.com', token: '']
        
        when:
        client = driver.createK8sClient([:])
        then:
        1 * driver.configDiscover(null) >> new ClientConfig(server: 'http://k8s.com:8000', token: 'xyz')
        client.config.server == 'http://k8s.com:8000'
        client.config.token == 'xyz'
        client.config.namespace == 'default'
        client.config.serviceAccount == null

        when:
        cfg = [k8s: [client:CLIENT_CFG, namespace: 'my-namespace', serviceAccount: 'my-account']]
        client = driver.createK8sClient(cfg)
        then:
        1 * driver.configCreate(CLIENT_CFG) >> { ClientConfig.fromMap(CLIENT_CFG) }
        client.config.server == 'foo.com'
        client.config.namespace == 'my-namespace'
        client.config.serviceAccount == 'my-account'

        when:
        cfg = [k8s: [context: 'foo']]
        client = driver.createK8sClient(cfg)
        then:
        1 * driver.configDiscover('foo') >> new ClientConfig(server: 'http://server.com', namespace: 'my-namespace')
        client.config.server == 'http://server.com'
        client.config.namespace == 'my-namespace'

    }

    def 'should create launcher spec' () {

        given:
        def driver = Spy(K8sDriverLauncher)

        when:
        driver.runName = 'foo-boo'
        driver.userDir = '/the/user/dir'
        driver.workDir = '/the/work/dir'
        driver.projectDir = '/the/project/dir'
        driver.configMounts['cfg-2'] = '/mnt/path/2'
        driver.client = new K8sClient(new ClientConfig(namespace: 'foo', serviceAccount: 'bar'))

        def spec = driver.makeLauncherSpec()
        then:
        driver.getImageName() >> 'the-image'
        driver.getVolumeClaims() >> new VolumeClaims( vol1: [mountPath: '/mnt/path/1'] )
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
                                         [name:'vol-1', mountPath:'/mnt/path/1'],
                                         [name:'vol-2', mountPath:'/mnt/path/2']]]
                                ],
                        serviceAccountName:'bar',
                        volumes:[[name:'vol-1', persistentVolumeClaim:[claimName:'vol1']],
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
        def SCM_CONFIG = [hello: 'world']

        def EXPECTED = [:]
        EXPECTED['init.sh'] == ''
        when:
        driver.userDir = '/launch/dir'
        driver.config = NXF_CONFIG
        driver.scm = SCM_CONFIG
        driver.cmd = new CmdKubeRun(paramsFile: params.toString())

        driver.createK8sConfigMap()
        then:
        1 * driver.makeConfigMapName(_ as Map) >> 'nf-config-123'
        1 * driver.tryCreateConfigMap('nf-config-123', _ as Map) >> {  name, cfg ->
            assert cfg.'init.sh' == "mkdir -p '/launch/dir'; if [ -d '/launch/dir' ]; then cd '/launch/dir'; else echo 'Cannot create nextflow userDir: /launch/dir'; exit 1; fi; [ -f /etc/nextflow/scm ] && ln -s /etc/nextflow/scm \$NXF_HOME/scm; [ -f /etc/nextflow/nextflow.config ] && cp /etc/nextflow/nextflow.config \$PWD/nextflow.config; echo cd \\\"\$PWD\\\" > /root/.profile; "
            assert cfg.'nextflow.config' == "foo = 'bar'\n"
            assert cfg.'scm' == "hello = 'world'\n"
            assert cfg.'params.json' == 'bla-bla'
            return null
        }

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
        CFG_WITH_MOUNTS.k8s.volumeClaims = [ pvc: [mountPath:'/foo'] ]

        when:
        driver.cmd = new CmdKubeRun()
        config = driver.makeConfig(NAME)
        then:
        1 *  driver.loadConfig(NAME) >> CFG_EMPTY
        config.process.executor == 'k8s'
        config.k8s.autoMountHostPaths == false

        when:
        driver.cmd = new CmdKubeRun()
        config = driver.makeConfig(NAME)
        then:
        1 *  driver.loadConfig(NAME) >> CFG_WITH_MOUNTS
        config.process.executor == 'k8s'
        config.k8s.autoMountHostPaths == false
        config.k8s.volumeClaims.pvc.mountPath == '/foo'
        config.k8s.volumeClaims.size() == 1

        when:
        driver.cmd = new CmdKubeRun(volMounts: ['pvc-1:/this','pvc-2:/that'] )
        config = driver.makeConfig(NAME)
        then:
        1 *  driver.loadConfig(NAME) >> CFG_EMPTY
        config.process.executor == 'k8s'
        config.k8s.autoMountHostPaths == false
        config.k8s.volumeClaims.'pvc-1'.mountPath == '/this'
        config.k8s.volumeClaims.'pvc-2'.mountPath == '/that'
        config.k8s.volumeClaims.size() == 2


        when:
        driver.cmd = new CmdKubeRun(volMounts: ['xyz:/this'] )
        config = driver.makeConfig(NAME)
        then:
        1 *  driver.loadConfig(NAME) >> CFG_WITH_MOUNTS
        config.process.executor == 'k8s'
        config.k8s.autoMountHostPaths == false
        config.k8s.volumeClaims.'xyz'.mountPath == '/this'
        config.k8s.volumeClaims.'pvc'.mountPath == '/foo'
        config.k8s.volumeClaims.size() == 2
        config.k8s.volumeClaims.iterator().next().key == 'xyz'  // the command line entry should be the first

    }
}
