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

package nextflow.container

import java.nio.file.Paths

import nextflow.util.MemoryUnit
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DockerBuilderTest extends Specification {


    def 'test docker mounts'() {

        given:
        def builder = [:] as DockerBuilder
        def files =  [Paths.get('/folder/data'),  Paths.get('/folder/db'), Paths.get('/folder/db') ]
        def real = [ Paths.get('/user/yo/nextflow/bin'), Paths.get('/user/yo/nextflow/work'), Paths.get('/db/pdb/local/data') ]
        def quotes =  [ Paths.get('/folder with blanks/A'), Paths.get('/folder with blanks/B') ]

        expect:
        builder.makeVolumes([]).toString() == '-v "$PWD":"$PWD"'
        builder.makeVolumes(files).toString() == '-v /folder:/folder -v "$PWD":"$PWD"'
        builder.makeVolumes(real).toString()  == '-v /user/yo/nextflow:/user/yo/nextflow -v /db/pdb/local/data:/db/pdb/local/data -v "$PWD":"$PWD"'
        builder.makeVolumes(quotes).toString() == '-v /folder\\ with\\ blanks:/folder\\ with\\ blanks -v "$PWD":"$PWD"'

    }


    def 'test docker env'() {

        given:
        def builder = [:] as DockerBuilder

        expect:
        builder.makeEnv('X=1').toString() == '-e "X=1"'
        builder.makeEnv([VAR_X:1, VAR_Y: 2]).toString() == '-e "VAR_X=1" -e "VAR_Y=2"'
        builder.makeEnv('BAR').toString() == '${BAR:+-e "BAR=$BAR"}'
    }

    def 'test docker create command line'() {

        setup:
        def env = [FOO: 1, BAR: 'hello world']
        def db_file = Paths.get('/home/db')

        expect:
        new DockerBuilder('fedora')
                .build()
                .runCommand == 'docker run -i -v "$PWD":"$PWD" -w "$PWD" fedora'

        new DockerBuilder('fedora')
                .addEnv(env)
                .build()
                .runCommand == 'docker run -i -e "FOO=1" -e "BAR=hello world" -v "$PWD":"$PWD" -w "$PWD" fedora'

        new DockerBuilder('ubuntu')
                .params(temp:'/hola')
                .build()
                .runCommand == 'docker run -i -v /hola:/tmp -v "$PWD":"$PWD" -w "$PWD" ubuntu'

        new DockerBuilder('busybox')
                .params(sudo: true)
                .build()
                .runCommand == 'sudo docker run -i -v "$PWD":"$PWD" -w "$PWD" busybox'

        new DockerBuilder('busybox')
                .params(entry: '/bin/bash')
                .build()
                .runCommand == 'docker run -i -v "$PWD":"$PWD" -w "$PWD" --entrypoint /bin/bash busybox'

        new DockerBuilder('busybox')
                .params(runOptions: '-x --zeta')
                .build()
                .runCommand == 'docker run -i -v "$PWD":"$PWD" -w "$PWD" -x --zeta busybox'

        new DockerBuilder('busybox')
                .params(userEmulation:true)
                .build()
                .runCommand == 'docker run -i -u $(id -u) -e "HOME=${HOME}" -v /etc/passwd:/etc/passwd:ro -v /etc/shadow:/etc/shadow:ro -v /etc/group:/etc/group:ro -v $HOME:$HOME -v "$PWD":"$PWD" -w "$PWD" busybox'

        new DockerBuilder('busybox')
                .setName('hola')
                .build()
                .runCommand == 'docker run -i -v "$PWD":"$PWD" -w "$PWD" --name hola busybox'

        new DockerBuilder('busybox')
                .params(engineOptions: '--tlsverify --tlscert="/path/to/my/cert"')
                .build()
                .runCommand == 'docker --tlsverify --tlscert="/path/to/my/cert" run -i -v "$PWD":"$PWD" -w "$PWD" busybox'

        new DockerBuilder('fedora')
                .addEnv(env)
                .addMount(db_file)
                .addMount(db_file)  // <-- add twice the same to prove that the final string won't contain duplicates
                .build()
                .runCommand == 'docker run -i -e "FOO=1" -e "BAR=hello world" -v /home/db:/home/db -v "$PWD":"$PWD" -w "$PWD" fedora'

        new DockerBuilder('fedora')
                .params(readOnlyInputs: true)
                .addMount(db_file)
                .build()
                .runCommand == 'docker run -i -v /home/db:/home/db:ro -v "$PWD":"$PWD" -w "$PWD" fedora'


        new DockerBuilder('fedora')
                .params(mountFlags: 'Z')
                .addMount(db_file)
                .build()
                .runCommand == 'docker run -i -v /home/db:/home/db:Z -v "$PWD":"$PWD":Z -w "$PWD" fedora'
    }

    def 'test memory and cpuset' () {

        expect:
        new DockerBuilder('fedora')
                .setCpus('1,2')
                .build()
                .runCommand == 'docker run -i --cpuset-cpus 1,2 -v "$PWD":"$PWD" -w "$PWD" fedora'

        new DockerBuilder('fedora')
                .params(legacy:true)
                .setCpus('1,2')
                .build()
                .runCommand == 'docker run -i --cpuset 1,2 -v "$PWD":"$PWD" -w "$PWD" fedora'

        new DockerBuilder('fedora')
                .setMemory('10g')
                .build()
                .runCommand == 'docker run -i --memory 10g -v "$PWD":"$PWD" -w "$PWD" fedora'

        new DockerBuilder('fedora')
                .setMemory(new MemoryUnit('100M'))
                .build()
                .runCommand == 'docker run -i --memory 100m -v "$PWD":"$PWD" -w "$PWD" fedora'

        new DockerBuilder('fedora')
                .setCpus('1-3')
                .setMemory(new MemoryUnit('100M'))
                .build()
                .runCommand == 'docker run -i --cpuset-cpus 1-3 --memory 100m -v "$PWD":"$PWD" -w "$PWD" fedora'

    }

    def 'test add mount'() {

        when:
        def docker = new DockerBuilder('fedora')
        docker.addMount(Paths.get('hello'))
        then:
        docker.mounts.size() == 1
        docker.mounts.contains(Paths.get('hello'))

        when:
        docker.addMount(null)
        then:
        docker.mounts.size() == 1

    }

    def 'test get commands'() {

        when:
        def docker =  new DockerBuilder('busybox').setName('c1').build()
        then:
        docker.runCommand == 'docker run -i -v "$PWD":"$PWD" -w "$PWD" --name c1 busybox'
        docker.removeCommand == 'docker rm c1'
        docker.killCommand == 'docker kill c1'

        when:
        docker =  new DockerBuilder('busybox').setName('c2').params(sudo: true, remove: true).build()
        then:
        docker.runCommand == 'sudo docker run -i -v "$PWD":"$PWD" -w "$PWD" --name c2 busybox'
        docker.removeCommand == 'sudo docker rm c2'
        docker.killCommand == 'sudo docker kill c2'


        when:
        docker =  new DockerBuilder('busybox').setName('c3').params(remove: true).build()
        then:
        docker.runCommand == 'docker run -i -v "$PWD":"$PWD" -w "$PWD" --name c3 busybox'
        docker.removeCommand == 'docker rm c3'
        docker.killCommand == 'docker kill c3'

        when:
        docker =  new DockerBuilder('busybox').setName('c4').params(kill: 'SIGKILL').build()
        then:
        docker.killCommand == 'docker kill -s SIGKILL c4'

        when:
        docker =  new DockerBuilder('busybox').setName('c5').params(kill: false,remove: false).build()
        then:
        docker.killCommand == null
        docker.removeCommand == null

    }

    def 'should add docker run to shell script' () {

        when:
        def script = '''
            #!/bin/bash
            FOO=bar
            busybox --foo --bar
            do_this
            do_that
            '''
        def tokens = ContainerScriptTokens.parse(script)
        def docker = new DockerBuilder('busybox').addEnv(tokens.variables)
        docker.build()

        then:
        docker.addContainerRunCommand(tokens) == '''
            #!/bin/bash
            FOO=bar
            docker run -i -e "FOO=bar" -v "$PWD":"$PWD" -w "$PWD" busybox --foo --bar
            do_this
            do_that
            '''
                .stripIndent().leftTrim()

        when:
        tokens = ContainerScriptTokens.parse('#!/bin/bash\nbusybox')
        docker = new DockerBuilder('busybox')
        docker.build()
        then:
        docker.addContainerRunCommand(tokens) == '''
            #!/bin/bash
            docker run -i -v "$PWD":"$PWD" -w "$PWD" busybox
            '''
                .stripIndent().leftTrim()


    }


    def 'should get run command line' () {

        when:
        def cli = new DockerBuilder('ubuntu:14').build().getRunCommand()
        then:
        cli ==  'docker run -i -v "$PWD":"$PWD" -w "$PWD" ubuntu:14'

        when:
        cli = new DockerBuilder('ubuntu:14').build().getRunCommand('bwa --this --that file.fasta')
        then:
        cli ==  'docker run -i -v "$PWD":"$PWD" -w "$PWD" ubuntu:14 bwa --this --that file.fasta'

        when:
        cli = new DockerBuilder('ubuntu:14').params(entry:'/bin/bash').build().getRunCommand('bwa --this --that file.fasta')
        then:
        cli ==  'docker run -i -v "$PWD":"$PWD" -w "$PWD" --entrypoint /bin/bash ubuntu:14 -c "bwa --this --that file.fasta"'

    }

    def 'should return mount flags'() {

        given:
        def builder = new DockerBuilder().params(mountFlags: flags)

        expect:
        builder.mountFlags(readOnly) == expected

        where:
        readOnly | flags    | expected
        false    | null     | ''
        false    | "Z"      | ':Z'
        false    | "z,Z "   | ':z,Z'
        true     | null     | ':ro'
        true     | ''       | ':ro'
        true     | 'Z'      | ':ro,Z'
    }


}
