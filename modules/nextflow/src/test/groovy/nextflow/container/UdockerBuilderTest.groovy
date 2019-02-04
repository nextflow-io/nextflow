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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class UdockerBuilderTest extends Specification {


    def 'test udocker env'() {

        given:
        def builder = new UdockerBuilder('x')

        expect:
        builder.makeEnv('X=1').toString() == '-e "X=1"'
        builder.makeEnv([VAR_X:1, VAR_Y: 2]).toString() == '-e "VAR_X=1" -e "VAR_Y=2"'
    }


    def 'test udocker mounts'() {

        given:
        def builder = new UdockerBuilder('x')
        def files =  [Paths.get('/folder/data'),  Paths.get('/folder/db'), Paths.get('/folder/db') ]
        def real = [ Paths.get('/user/yo/nextflow/bin'), Paths.get('/user/yo/nextflow/work'), Paths.get('/db/pdb/local/data') ]
        def quotes =  [ Paths.get('/folder with blanks/A'), Paths.get('/folder with blanks/B') ]

        expect:
        builder.makeVolumes([]).toString() == '-v "$PWD":"$PWD"'
        builder.makeVolumes(files).toString() == '-v /folder:/folder -v "$PWD":"$PWD"'
        builder.makeVolumes(real).toString()  == '-v /user/yo/nextflow:/user/yo/nextflow -v /db/pdb/local/data:/db/pdb/local/data -v "$PWD":"$PWD"'
        builder.makeVolumes(quotes).toString() == '-v /folder\\ with\\ blanks:/folder\\ with\\ blanks -v "$PWD":"$PWD"'

    }

    def 'should get run cmd line' () {

        given:
        def env = [FOO:1, BAR:'hello world']
        def db_file = Paths.get('/home/db')

        expect:
        new UdockerBuilder('fedora')
                .build()
                .@runCommand == 'udocker.py run --rm -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "fedora:latest")'

        new UdockerBuilder('fedora')
                .addEnv(env)
                .build()
                .@runCommand == 'udocker.py run --rm -e "FOO=1" -e "BAR=hello world" -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "fedora:latest")'

        new UdockerBuilder('fedora')
                .setCpus('1,2')
                .build()
                .@runCommand == 'udocker.py run --rm --cpuset-cpus=1,2 -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "fedora:latest")'

        new UdockerBuilder('fedora')
                .addMount(db_file)
                .addEnv(env)
                .build()
                .@runCommand == 'udocker.py run --rm -e "FOO=1" -e "BAR=hello world" -v /home/db:/home/db -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "fedora:latest")'

        new UdockerBuilder('busybox')
                .params(remove: false)
                .build()
                .@runCommand == 'udocker.py run -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "busybox:latest")'

        new UdockerBuilder('busybox')
                .params(runOptions: '-x --zeta')
                .build()
                .@runCommand == 'udocker.py run --rm -v "$PWD":"$PWD" -w "$PWD" --bindhome -x --zeta $(udocker.py create "busybox:latest")'

        new UdockerBuilder('busybox')
                .params(entry: '/bin/blah')
                .build()
                .@runCommand == 'udocker.py run --rm -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "busybox:latest")'

    }

    def 'should append the run command line' () {

        given:
        def builder = new UdockerBuilder('ubuntu:latest')

        when:
        def result = builder.build().getRunCommand()
        then:
        result == '''
            ((udocker.py images | egrep -o "^ubuntu:latest\\s") || udocker.py pull "ubuntu:latest")>/dev/null
            [[ $? != 0 ]] && echo "Udocker failed while pulling container \\`ubuntu:latest\\`" >&2 && exit 1
            udocker.py run --rm -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "ubuntu:latest")
            '''
            .stripIndent().trim()

        builder.getRemoveCommand() == null
        builder.getKillCommand() == null
    }

    def 'should append the run command line with launcher' () {

        when:
        def builder = new UdockerBuilder('ubuntu:latest')
        def result = builder.build().getRunCommand('bwa --this --that')
        then:
        result == '''
            ((udocker.py images | egrep -o "^ubuntu:latest\\s") || udocker.py pull "ubuntu:latest")>/dev/null
            [[ $? != 0 ]] && echo "Udocker failed while pulling container \\`ubuntu:latest\\`" >&2 && exit 1
            udocker.py run --rm -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "ubuntu:latest") bwa --this --that
            '''
                .stripIndent().trim()

        builder.getRemoveCommand() == null
        builder.getKillCommand() == null


        when:
        builder = new UdockerBuilder('ubuntu:latest').params(entry:'/bin/bash')
        result = builder.build().getRunCommand('bwa --this --that')
        then:
        result == '''
            ((udocker.py images | egrep -o "^ubuntu:latest\\s") || udocker.py pull "ubuntu:latest")>/dev/null
            [[ $? != 0 ]] && echo "Udocker failed while pulling container \\`ubuntu:latest\\`" >&2 && exit 1
            udocker.py run --rm -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "ubuntu:latest") /bin/bash -c "bwa --this --that"
            '''
                .stripIndent().trim()

        builder.getRemoveCommand() == null
        builder.getKillCommand() == null
    }

}
