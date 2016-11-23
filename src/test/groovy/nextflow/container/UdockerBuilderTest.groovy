/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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
        builder.makeEnv( Paths.get('/some/file.env') ).toString() == '-e "BASH_ENV=/some/file.env"'
        builder.makeEnv( new File('/some/file.env') ).toString() == '-e "BASH_ENV=/some/file.env"'
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
        def envFile = Paths.get('/data/env.file')
        def db_file = Paths.get('/home/db')

        expect:
        new UdockerBuilder('fedora')
                .build()
                .runCommand == 'udocker.py run --rm -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "fedora:latest") /bin/bash'

        new UdockerBuilder('fedora')
                .addEnv(envFile)
                .build()
                .runCommand == 'udocker.py run --rm -e "BASH_ENV=/data/env.file" -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "fedora:latest") /bin/bash'

        new UdockerBuilder('fedora')
                .setCpus('1,2')
                .build()
                .runCommand == 'udocker.py run --rm --cpuset-cpus=1,2 -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "fedora:latest") /bin/bash'

        new UdockerBuilder('fedora')
                .addMount(db_file)
                .addEnv(envFile)
                .build()
                .runCommand == 'udocker.py run --rm -e "BASH_ENV=/data/env.file" -v /home/db:/home/db -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "fedora:latest") /bin/bash'

        new UdockerBuilder('busybox')
                .params(remove: false)
                .build()
                .runCommand == 'udocker.py run -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "busybox:latest") /bin/bash'

        new UdockerBuilder('busybox')
                .params(runOptions: '-x --zeta')
                .build()
                .runCommand == 'udocker.py run --rm -v "$PWD":"$PWD" -w "$PWD" --bindhome -x --zeta $(udocker.py create "busybox:latest") /bin/bash'

        new UdockerBuilder('busybox')
                .params(entry: '/bin/blah')
                .build()
                .runCommand == 'udocker.py run --rm -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "busybox:latest") /bin/blah'

    }

    def 'should append the run command line' () {

        given:
        def script = new StringBuilder()

        when:
        def builder = new UdockerBuilder('ubuntu:latest')
        builder.build()
        builder.appendRunCommand(script)
        then:
        script.toString() == '''
            ((udocker.py images | egrep -o "^ubuntu:latest\\s") || udocker.py pull "ubuntu:latest")>/dev/null
            [[ $? != 0 ]] && echo "Udocker failed while pulling container \\`ubuntu:latest\\`" >&2 && exit 1
            udocker.py run --rm -v "$PWD":"$PWD" -w "$PWD" --bindhome $(udocker.py create "ubuntu:latest") /bin/bash
            '''
            .stripIndent().trim()

        builder.getRemoveCommand() == null
        builder.getKillCommand() == null
    }

}
