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

class AppleContainerBuilderTest extends Specification {

    def 'should build the basic run command'() {
        expect:
        new AppleContainerBuilder('alpine')
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" alpine'
    }

    def 'should include env vars with -e'() {
        expect:
        new AppleContainerBuilder('alpine')
                .addEnv([FOO: 1, BAR: 'hello world'])
                .build()
                .runCommand == 'container run -i -e "FOO=1" -e "BAR=hello world" -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" alpine'
    }

    def 'should include mounts'() {
        given:
        def db = Paths.get('/home/db')
        expect:
        new AppleContainerBuilder('alpine')
                .addMount(db)
                .build()
                .runCommand == 'container run -i -v /home/db:/home/db -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" alpine'
    }

    def 'should mount temp dir'() {
        expect:
        new AppleContainerBuilder('alpine', new AppleContainerConfig(temp: '/hola'))
                .build()
                .runCommand == 'container run -i -v /hola:/tmp -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" alpine'
    }

    def 'should add -t when tty enabled'() {
        expect:
        new AppleContainerBuilder('alpine', new AppleContainerConfig(tty: true))
                .build()
                .runCommand == 'container run -i -t -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" alpine'
    }

    def 'should pass cpus as --cpus'() {
        expect:
        new AppleContainerBuilder('alpine')
                .setCpus(4)
                .build()
                .runCommand == 'container run -i --cpus 4 -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" alpine'
    }

    def 'should pass memory as -m'() {
        expect:
        new AppleContainerBuilder('alpine')
                .setMemory('2G')
                .build()
                .runCommand == 'container run -i -m 2G -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" alpine'

        new AppleContainerBuilder('alpine')
                .setMemory(new MemoryUnit('100M'))
                .build()
                .runCommand == 'container run -i -m 100m -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" alpine'
    }

    def 'should pass platform as --platform'() {
        expect:
        new AppleContainerBuilder('alpine')
                .setPlatform('linux/arm64')
                .build()
                .runCommand == 'container run -i --platform linux/arm64 -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" alpine'
    }

    def 'should emit entrypoint via params'() {
        expect:
        new AppleContainerBuilder('alpine')
                .params(entry: '/bin/bash')
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --entrypoint /bin/bash alpine'
    }

    def 'should emit cap-add via params'() {
        expect:
        new AppleContainerBuilder('alpine')
                .params(capAdd: 'SYS_ADMIN')
                .build()
                .runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --cap-add SYS_ADMIN alpine'
    }

    def 'should set name and produce remove/kill commands'() {
        when:
        def b = new AppleContainerBuilder('alpine').setName('c1').build()
        then:
        b.runCommand == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --name c1 alpine'
        b.removeCommand == 'container rm c1'
        b.killCommand == 'container stop c1'
    }

    def 'should disable remove and kill when configured'() {
        when:
        def config = new AppleContainerConfig(remove: false)
        def b = new AppleContainerBuilder('alpine', config).setName('c2').params(kill: false).build()
        then:
        b.removeCommand == null
        b.killCommand == null
    }

    def 'should emit kill with signal string'() {
        when:
        def b = new AppleContainerBuilder('alpine').setName('c3').params(kill: 'SIGKILL').build()
        then:
        b.killCommand == 'container kill -s SIGKILL c3'
    }

    def 'should append runOptions and engineOptions'() {
        expect:
        new AppleContainerBuilder('alpine', new AppleContainerConfig(runOptions: '--ssh', engineOptions: '--debug'))
                .build()
                .runCommand == 'container --debug run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --ssh alpine'
    }

    def 'should append launcher after run command'() {
        when:
        def cli = new AppleContainerBuilder('alpine').build().getRunCommand('bwa --this file.fa')
        then:
        cli == 'container run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" alpine bwa --this file.fa'
    }
}
