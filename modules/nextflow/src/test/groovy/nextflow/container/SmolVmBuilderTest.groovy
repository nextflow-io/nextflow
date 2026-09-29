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

package nextflow.container

import java.nio.file.Paths

import nextflow.util.MemoryUnit
import spock.lang.Specification

class SmolVmBuilderTest extends Specification {

    def 'should build the basic run command'() {
        expect:
        new SmolVmBuilder('alpine')
                .build()
                .runCommand == 'smolvm machine run -i --net -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --image alpine'
    }

    def 'should include env vars with -e'() {
        expect:
        new SmolVmBuilder('alpine')
                .addEnv([FOO: 1, BAR: 'hello world'])
                .build()
                .runCommand == 'smolvm machine run -i --net -e "FOO=1" -e "BAR=hello world" -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --image alpine'
    }

    def 'should include mounts'() {
        given:
        def db = Paths.get('/home/db')
        expect:
        new SmolVmBuilder('alpine')
                .addMount(db)
                .build()
                .runCommand == 'smolvm machine run -i --net -v /home/db:/home/db -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --image alpine'
    }

    def 'should mount inputs read-only when configured'() {
        given:
        def db = Paths.get('/home/db')
        expect:
        new SmolVmBuilder('alpine', new SmolVmConfig(writableInputMounts: false))
                .addMount(db)
                .build()
                .runCommand == 'smolvm machine run -i --net -v /home/db:/home/db:ro -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --image alpine'
    }

    def 'should pass cpus as --cpus'() {
        expect:
        new SmolVmBuilder('alpine')
                .setCpus(4)
                .build()
                .runCommand == 'smolvm machine run -i --net --cpus 4 -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --image alpine'
    }

    def 'should pass memory as --mem in MiB'() {
        expect:
        new SmolVmBuilder('alpine')
                .setMemory(new MemoryUnit('2G'))
                .build()
                .runCommand == 'smolvm machine run -i --net --mem 2048 -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --image alpine'

        new SmolVmBuilder('alpine')
                .setMemory('512m')
                .build()
                .runCommand == 'smolvm machine run -i --net --mem 512 -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --image alpine'
    }

    def 'should pass platform as --oci-platform'() {
        expect:
        new SmolVmBuilder('alpine')
                .setPlatform('linux/arm64')
                .build()
                .runCommand == 'smolvm machine run -i --net --oci-platform linux/arm64 -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --image alpine'
    }

    def 'should add -t when tty enabled'() {
        expect:
        new SmolVmBuilder('alpine', new SmolVmConfig(tty: true))
                .build()
                .runCommand == 'smolvm machine run -i -t --net -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --image alpine'
    }

    def 'should omit --net when network disabled'() {
        expect:
        new SmolVmBuilder('alpine', new SmolVmConfig(network: false))
                .build()
                .runCommand == 'smolvm machine run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --image alpine'
    }

    def 'should append runOptions and engineOptions'() {
        expect:
        new SmolVmBuilder('alpine', new SmolVmConfig(runOptions: '--ssh-agent', engineOptions: '--verbose'))
                .build()
                .runCommand == 'smolvm --verbose machine run -i --net -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --ssh-agent --image alpine'
    }

    def 'should mount temp dir'() {
        expect:
        new SmolVmBuilder('alpine', new SmolVmConfig(temp: '/hola'))
                .build()
                .runCommand == 'smolvm machine run -i --net -v /hola:/tmp -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --image alpine'
    }

    def 'should append launcher after run command'() {
        when:
        def cli = new SmolVmBuilder('alpine').build().getRunCommand('bwa --this file.fa')
        then:
        cli == 'smolvm machine run -i --net -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --image alpine bwa --this file.fa'
    }

    def 'should not provide remove or kill commands (ephemeral VM)'() {
        when:
        def b = new SmolVmBuilder('alpine').build()
        then:
        b.removeCommand == null
        // base default kill command (exit signal), not an smolvm command
        b.killCommand != null
        !b.killCommand.contains('smolvm')
    }
}
