/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
 * @author tbugfinder <github@online.ms>
 */
class PodmanBuilderTest extends Specification {


    def 'test podman mounts'() {

        given:
        def builder = [:] as PodmanBuilder
        def files =  [Paths.get('/folder/data'),  Paths.get('/folder/db'), Paths.get('/folder/db') ]
        def real = [ Paths.get('/user/yo/nextflow/bin'), Paths.get('/user/yo/nextflow/work'), Paths.get('/db/pdb/local/data') ]
        def quotes =  [ Paths.get('/folder with blanks/A'), Paths.get('/folder with blanks/B') ]

        expect:
        builder.makeVolumes([]).toString() == '-v "$PWD":"$PWD"'
        builder.makeVolumes(files).toString() == '-v /folder:/folder -v "$PWD":"$PWD"'
        builder.makeVolumes(real).toString()  == '-v /user/yo/nextflow:/user/yo/nextflow -v /db/pdb/local/data:/db/pdb/local/data -v "$PWD":"$PWD"'
        builder.makeVolumes(quotes).toString() == '-v /folder\\ with\\ blanks:/folder\\ with\\ blanks -v "$PWD":"$PWD"'

    }


    def 'test podman env'() {

        given:
        def builder = [:] as PodmanBuilder

        expect:
        builder.makeEnv('X=1').toString() == '-e "X=1"'
        builder.makeEnv([VAR_X:1, VAR_Y: 2]).toString() == '-e "VAR_X=1" -e "VAR_Y=2"'
        builder.makeEnv('BAR').toString() == '${BAR:+-e "BAR=$BAR"}'
    }

    def 'test podman create command line'() {

        setup:
        def env = [FOO: 1, BAR: 'hello world']
        def db_file = Paths.get('/home/db')

        expect:
        new PodmanBuilder('fedora')
                .build()
                .runCommand == 'podman run -i -v "$PWD":"$PWD" -w "$PWD" fedora'

        new PodmanBuilder('fedora')
                .addEnv(env)
                .build()
                .runCommand == 'podman run -i -e "FOO=1" -e "BAR=hello world" -v "$PWD":"$PWD" -w "$PWD" fedora'

        new PodmanBuilder('ubuntu')
                .params(temp:'/hola')
                .build()
                .runCommand == 'podman run -i -v /hola:/tmp -v "$PWD":"$PWD" -w "$PWD" ubuntu'

        new PodmanBuilder('busybox')
                .params(entry: '/bin/bash')
                .build()
                .runCommand == 'podman run -i -v "$PWD":"$PWD" -w "$PWD" --entrypoint /bin/bash busybox'

        new PodmanBuilder('busybox')
                .params(runOptions: '-x --zeta')
                .build()
                .runCommand == 'podman run -i -v "$PWD":"$PWD" -w "$PWD" -x --zeta busybox'

        new PodmanBuilder('busybox')
                .setName('hola')
                .build()
                .runCommand == 'podman run -i -v "$PWD":"$PWD" -w "$PWD" --name hola busybox'

        new PodmanBuilder('busybox')
                .params(engineOptions: '--tls-verify=false --cert-dir "/path/to/my/cert-dir"')
                .build()
                .runCommand == 'podman --tls-verify=false --cert-dir "/path/to/my/cert-dir" run -i -v "$PWD":"$PWD" -w "$PWD" busybox'

        new PodmanBuilder('fedora')
                .addEnv(env)
                .addMount(db_file)
                .addMount(db_file)  // <-- add twice the same to prove that the final string won't contain duplicates
                .build()
                .runCommand == 'podman run -i -e "FOO=1" -e "BAR=hello world" -v /home/db:/home/db -v "$PWD":"$PWD" -w "$PWD" fedora'

        new PodmanBuilder('fedora')
                .params(readOnlyInputs: true)
                .addMount(db_file)
                .build()
                .runCommand == 'podman run -i -v /home/db:/home/db:ro -v "$PWD":"$PWD" -w "$PWD" fedora'


        new PodmanBuilder('fedora')
                .params(mountFlags: 'Z')
                .addMount(db_file)
                .build()
                .runCommand == 'podman run -i -v /home/db:/home/db:Z -v "$PWD":"$PWD":Z -w "$PWD" fedora'
    }

    def 'test add mount'() {

        when:
        def podman = new PodmanBuilder('fedora')
        podman.addMount(Paths.get('hello'))
        then:
        podman.mounts.size() == 1
        podman.mounts.contains(Paths.get('hello'))

        when:
        podman.addMount(null)
        then:
        podman.mounts.size() == 1

    }

    def 'test get commands'() {

        when:
        def podman = new PodmanBuilder('busybox').setName('c1').build()
        then:
        podman.runCommand == 'podman run -i -v "$PWD":"$PWD" -w "$PWD" --name c1 busybox'
        podman.removeCommand == 'podman rm c1'
        podman.killCommand == 'podman kill c1'

        when:
        podman = new PodmanBuilder('busybox').setName('c3').params(remove: true).build()
        then:
        podman.runCommand == 'podman run -i -v "$PWD":"$PWD" -w "$PWD" --name c3 busybox'
        podman.removeCommand == 'podman rm c3'
        podman.killCommand == 'podman kill c3'

        when:
        podman = new PodmanBuilder('busybox').setName('c4').params(kill: 'SIGKILL').build()
        then:
        podman.killCommand == 'podman kill -s SIGKILL c4'

        when:
        podman = new PodmanBuilder('busybox').setName('c5').params(kill: false,remove: false).build()
        then:
        podman.killCommand == null
        podman.removeCommand == null

    }


    def 'should get run command line' () {

        when:
        def cli = new PodmanBuilder('ubuntu:14').build().getRunCommand()
        then:
        cli ==  'podman run -i -v "$PWD":"$PWD" -w "$PWD" ubuntu:14'

        when:
        cli = new PodmanBuilder('ubuntu:14').build().getRunCommand('bwa --this --that file.fasta')
        then:
        cli ==  'podman run -i -v "$PWD":"$PWD" -w "$PWD" ubuntu:14 bwa --this --that file.fasta'

        when:
        cli = new PodmanBuilder('ubuntu:14').params(entry:'/bin/bash').build().getRunCommand('bwa --this --that file.fasta')
        then:
        cli ==  'podman run -i -v "$PWD":"$PWD" -w "$PWD" --entrypoint /bin/bash ubuntu:14 -c "bwa --this --that file.fasta"'

    }

    def 'should return mount flags'() {

        given:
        def builder = new PodmanBuilder().params(mountFlags: flags)

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
