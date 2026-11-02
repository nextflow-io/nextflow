/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.executor.local

import java.nio.file.Path
import java.nio.file.Paths

import nextflow.SysEnv
import nextflow.container.DockerConfig
import nextflow.container.PodmanConfig
import nextflow.processor.TaskBean
import nextflow.util.MemoryUnit
import spock.lang.Specification

class DockerWrapperLauncherTest extends Specification {

    def cleanup() {
        SysEnv.pop()
    }

    def setup() {
        SysEnv.push([:])
    }

    private TaskBean bean(Map opts=[:]) {
        final result = new TaskBean(
            name: 'foo',
            workDir: Path.of('/work/ab/cdef'),
            containerImage: 'ubuntu:24.04',
            containerConfig: new DockerConfig(enabled: true),
            inputFiles: [:] )
        opts.each { k, v -> result.setProperty(k as String, v) }
        return result
    }

    def 'should create the default run command' () {
        given:
        def launcher = new DockerWrapperLauncher(bean(), 'nxf-1234')

        expect:
        launcher.name == 'nxf-1234'
        launcher.runCommand() == [
            'docker', 'run', '--rm',
            '-v', '/work/ab/cdef:/work/ab/cdef',
            '-w', '/work/ab/cdef',
            '--name', 'nxf-1234',
            'ubuntu:24.04',
            '/bin/bash', '-ue', '.command.run' ]
    }

    def 'should mount inputs and bin dirs' () {
        given:
        def task = bean(
            inputFiles: [
                'a.fq': Path.of('/data/reads/a.fq'),
                'b.fq': Path.of('/data/reads/b.fq'),
                'c.txt': Path.of('/data/c.txt'),
                'd.txt': Path.of('/work/xy/z123/d.txt') ],
            binDirs: [Path.of('/project/bin')] )
        def launcher = new DockerWrapperLauncher(task, 'nxf-1234')

        expect:
        launcher.volumes() == [
            '-v', '/data:/data',
            '-v', '/project/bin:/project/bin',
            '-v', '/work/ab/cdef:/work/ab/cdef',
            '-v', '/work/xy/z123:/work/xy/z123' ]
    }

    def 'should mount inputs read-only and the work dir writable' () {
        given:
        def config = new DockerConfig(enabled: true, writableInputMounts: false, mountFlags: 'Z')
        def task = bean(
            containerConfig: config,
            workDir: Path.of('/project/work/ab/cdef'),
            inputFiles: [
                'a.fq': Path.of('/project/data/a.fq'),
                'b.fq': Path.of('/project/b.fq'),
                'c.txt': Path.of('/project/work/ab/cdef/sub/c.txt') ] )
        def launcher = new DockerWrapperLauncher(task, 'nxf-1234')

        expect:
        launcher.volumes() == [
            '-v', '/project:/project:ro,Z',
            '-v', '/project/work/ab/cdef:/project/work/ab/cdef:Z' ]
    }

    def 'should map windows paths' () {
        given:
        def task = bean(
            workDir: Paths.get('C:\\Users\\me\\work\\ab\\cdef'),
            binDirs: [Paths.get('C:\\Users\\me\\My Project\\bin')] )
        def launcher = new DockerWrapperLauncher(task, 'nxf-1234')

        expect:
        launcher.runCommand() == [
            'docker', 'run', '--rm',
            '-v', 'C:\\Users\\me\\My Project\\bin:/c/Users/me/My Project/bin',
            '-v', 'C:\\Users\\me\\work\\ab\\cdef:/c/Users/me/work/ab/cdef',
            '-w', '/c/Users/me/work/ab/cdef',
            '--name', 'nxf-1234',
            'ubuntu:24.04',
            '/bin/bash', '-ue', '.command.run' ]
    }

    def 'should add resources, environment and user options' () {
        given:
        SysEnv.push(NXF_DEBUG: '1')
        def config = new DockerConfig(
            enabled: true,
            sudo: true,
            remove: false,
            engineOptions: '--context remote',
            runOptions: '--gpus all --label "x=hello world"',
            envWhitelist: 'FOO,BAR=1',
            temp: '/scratch/tmp' )
        def task = bean(
            containerConfig: config,
            containerCpus: 2,
            containerMemory: MemoryUnit.of('1 GB'),
            containerPlatform: 'linux/amd64',
            containerOptions: "--shm-size '2g'" )
        def launcher = new DockerWrapperLauncher(task, 'nxf-1234')

        expect:
        launcher.runCommand() == [
            'sudo', 'docker', '--context', 'remote', 'run',
            '--cpu-shares', '2048',
            '--memory', '1024m',
            '--platform', 'linux/amd64',
            '-e', 'NXF_DEBUG',
            '-e', 'FOO',
            '-e', 'BAR=1',
            '-v', '/scratch/tmp:/tmp',
            '-v', '/work/ab/cdef:/work/ab/cdef',
            '-w', '/work/ab/cdef',
            '--gpus', 'all', '--label', 'x=hello world',
            '--shm-size', '2g',
            '--name', 'nxf-1234',
            'ubuntu:24.04',
            '/bin/bash', '-ue', '.command.run' ]

        cleanup:
        SysEnv.pop()
    }

    def 'should use cpu limits and override the entrypoint' () {
        given:
        SysEnv.push(NXF_CONTAINER_ENTRYPOINT_OVERRIDE: 'true')
        def task = bean(
            containerConfig: new DockerConfig(enabled: true, cpuLimits: true, temp: 'auto'),
            containerCpus: 4 )
        def launcher = new DockerWrapperLauncher(task, 'nxf-1234')

        expect:
        launcher.runCommand() == [
            'docker', 'run', '--rm',
            '--cpus', '4',
            '-v', '/work/ab/cdef:/work/ab/cdef',
            '-w', '/work/ab/cdef',
            '--entrypoint', '/bin/bash',
            '--name', 'nxf-1234',
            'ubuntu:24.04',
            '-ue', '.command.run' ]

        cleanup:
        SysEnv.pop()
    }

    def 'should create the stop command' () {
        expect:
        new DockerWrapperLauncher(bean(containerConfig: CONFIG), 'nxf-1234').stopCommand() == EXPECTED

        where:
        CONFIG                                                          | EXPECTED
        new DockerConfig(enabled: true)                                 | ['docker', 'stop', 'nxf-1234']
        new DockerConfig(enabled: true, kill: 'SIGKILL')                | ['docker', 'kill', '-s', 'SIGKILL', 'nxf-1234']
        new DockerConfig(enabled: true, kill: false)                    | null
        new DockerConfig(enabled: true, sudo: true, engineOptions: '-H tcp://host:2375') | ['sudo', 'docker', '-H', 'tcp://host:2375', 'stop', 'nxf-1234']
    }

    def 'should reject missing image or non docker engine' () {
        when:
        new DockerWrapperLauncher(bean(containerImage: null), 'nxf-1234')
        then:
        thrown(IllegalArgumentException)

        when:
        new DockerWrapperLauncher(bean(containerConfig: new PodmanConfig(enabled: true)), 'nxf-1234')
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('requires the Docker engine')
    }

}
