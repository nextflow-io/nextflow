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

package nextflow.container

import java.nio.file.Paths

import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DockerBuilderTest extends Specification {


    def 'test docker mounts'() {

        given:
        def builder = Spy(DockerBuilder)
        def files =  [Paths.get('/folder/data'),  Paths.get('/folder/db'), Paths.get('/folder/db') ]
        def real = [ Paths.get('/user/yo/nextflow/bin'), Paths.get('/user/yo/nextflow/work'), Paths.get('/db/pdb/local/data') ]
        def quotes =  [ Paths.get('/folder with blanks/A'), Paths.get('/folder with blanks/B') ]

        expect:
        builder.makeVolumes([]).toString() == '-v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" '
        builder.makeVolumes(files).toString() == '-v /folder:/folder -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" '
        builder.makeVolumes(real).toString()  == '-v /user/yo/nextflow:/user/yo/nextflow -v /db/pdb/local/data:/db/pdb/local/data -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" '
        builder.makeVolumes(quotes).toString() == '-v /folder\\ with\\ blanks:/folder\\ with\\ blanks -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" '

    }

    @Unroll
    def 'test docker env'() {

        given:
        def builder = Spy(DockerBuilder)

        expect:
        builder.makeEnv(ENV).toString() == EXPECT

        where:
        ENV                 | EXPECT
        'X=1'               | '-e "X=1"'
        [VAR_X:1, VAR_Y: 2] | '-e "VAR_X=1" -e "VAR_Y=2"'
        'BAR'               | '-e "BAR"'
    }

    def 'test docker create command line'() {

        setup:
        def env = [FOO: 1, BAR: 'hello world']
        def db_file = Paths.get('/home/db')

        expect:
        new DockerBuilder('fedora')
                .build()
                .runCommand == 'docker run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'

        new DockerBuilder('fedora')
                .addEnv(env)
                .build()
                .runCommand == 'docker run -i -e "FOO=1" -e "BAR=hello world" -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'

        new DockerBuilder('ubuntu')
                .params(temp:'/hola')
                .build()
                .runCommand == 'docker run -i -v /hola:/tmp -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" ubuntu'

        new DockerBuilder('busybox')
                .params(sudo: true)
                .build()
                .runCommand == 'sudo docker run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" busybox'

        new DockerBuilder('busybox')
                .params(entry: '/bin/bash')
                .build()
                .runCommand == 'docker run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --entrypoint /bin/bash busybox'

        new DockerBuilder('busybox')
                .params(runOptions: '-x --zeta')
                .build()
                .runCommand == 'docker run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" -x --zeta busybox'

        new DockerBuilder('busybox')
                .setName('hola')
                .build()
                .runCommand == 'docker run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --name hola busybox'

        new DockerBuilder('busybox')
                .params(engineOptions: '--tlsverify --tlscert="/path/to/my/cert"')
                .build()
                .runCommand == 'docker --tlsverify --tlscert="/path/to/my/cert" run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" busybox'

        new DockerBuilder('fedora')
                .addEnv(env)
                .addMount(db_file)
                .addMount(db_file)  // <-- add twice the same to prove that the final string won't contain duplicates
                .build()
                .runCommand == 'docker run -i -e "FOO=1" -e "BAR=hello world" -v /home/db:/home/db -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'

        new DockerBuilder('fedora')
                .params(readOnlyInputs: true)
                .addMount(db_file)
                .build()
                .runCommand == 'docker run -i -v /home/db:/home/db:ro -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'


        new DockerBuilder('fedora')
                .params(mountFlags: 'Z')
                .addMount(db_file)
                .build()
                .runCommand == 'docker run -i -v /home/db:/home/db:Z -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR":Z -w "$NXF_TASK_WORKDIR" fedora'
    }

    def 'test memory and cpuset' () {

        expect:
        new DockerBuilder('fedora')
                .setCpus(2)
                .build()
                .runCommand == 'docker run -i --cpu-shares 2048 -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'

        new DockerBuilder('fedora')
                .setCpus(1)
                .build()
                .runCommand == 'docker run -i --cpu-shares 1024 -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'

        new DockerBuilder('fedora')
                .setCpus(8)
                .setCpuset('1,2')
                .build()
                .runCommand == 'docker run -i --cpu-shares 8192 --cpuset-cpus 1,2 -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'

        new DockerBuilder('fedora')
                .params(legacy: true)
                .setCpuset('1,2')
                .build()
                .runCommand == 'docker run -i --cpuset 1,2 -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'


        new DockerBuilder('fedora')
                .params(legacy: true)
                .setCpus(1)
                .build()
                .runCommand == 'docker run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'

        new DockerBuilder('fedora')
                .setMemory('10g')
                .build()
                .runCommand == 'docker run -i --memory 10g -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'

        new DockerBuilder('fedora')
                .setMemory(new MemoryUnit('100M'))
                .build()
                .runCommand == 'docker run -i --memory 100m -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'

        new DockerBuilder('fedora')
                .setCpuset('1-3')
                .setMemory(new MemoryUnit('100M'))
                .build()
                .runCommand == 'docker run -i --cpuset-cpus 1-3 --memory 100m -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" fedora'

        new DockerBuilder('fedora')
                .params(privileged: true)
                .build()
                .runCommand == 'docker run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --privileged fedora'

        new DockerBuilder('fedora')
                .params(device: '/dev/fuse')
                .params(capAdd: 'SYS_ADMIN')
                .build()
                .runCommand == 'docker run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --device /dev/fuse --cap-add SYS_ADMIN fedora'
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
        docker.runCommand == 'docker run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --name c1 busybox'
        docker.removeCommand == 'docker rm c1'
        docker.killCommand == 'docker stop c1'

        when:
        docker =  new DockerBuilder('busybox').setName('c2').params(sudo: true, remove: true).build()
        then:
        docker.runCommand == 'sudo docker run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --name c2 busybox'
        docker.removeCommand == 'sudo docker rm c2'
        docker.killCommand == 'sudo docker stop c2'


        when:
        docker =  new DockerBuilder('busybox').setName('c3').params(remove: true).build()
        then:
        docker.runCommand == 'docker run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --name c3 busybox'
        docker.removeCommand == 'docker rm c3'
        docker.killCommand == 'docker stop c3'

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


    def 'should get run command line' () {

        when:
        def cli = new DockerBuilder('ubuntu:14').build().getRunCommand()
        then:
        cli ==  'docker run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" ubuntu:14'

        when:
        cli = new DockerBuilder('ubuntu:14').build().getRunCommand('bwa --this --that file.fasta')
        then:
        cli ==  'docker run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" ubuntu:14 bwa --this --that file.fasta'

        when:
        cli = new DockerBuilder('ubuntu:14').params(entry:'/bin/bash').build().getRunCommand('bwa --this --that file.fasta')
        then:
        cli ==  'docker run -i -v "$NXF_TASK_WORKDIR":"$NXF_TASK_WORKDIR" -w "$NXF_TASK_WORKDIR" --entrypoint /bin/bash ubuntu:14 -c "bwa --this --that file.fasta"'

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
