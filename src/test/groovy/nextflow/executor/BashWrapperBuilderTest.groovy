/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.executor

import java.nio.file.Files
import java.nio.file.Paths

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BashWrapperBuilderTest extends Specification {

    def 'test changeToScratchDir' () {

        setup:
        def builder = [:] as BashWrapperBuilder

        expect:
        builder.changeToScratchDirectory() == null

        when:
        builder.scratch = true
        then:
        builder.changeToScratchDirectory() == 'NXF_SCRATCH=${TMPDIR:-`mktemp -d`} && cd $NXF_SCRATCH'

        when:
        builder.scratch = '$SOME_DIR'
        then:
        builder.changeToScratchDirectory() == 'NXF_SCRATCH=${SOME_DIR:-`mktemp -d`} && cd $NXF_SCRATCH'

        when:
        builder.scratch = '/my/temp'
        then:
        builder.changeToScratchDirectory() == 'NXF_SCRATCH=$(mktemp -d -p /my/temp) && cd $NXF_SCRATCH'

        when:
        builder.scratch = '/my/temp'
        then:
        builder.changeToScratchDirectory() == 'NXF_SCRATCH=$(mktemp -d -p /my/temp) && cd $NXF_SCRATCH'

    }


    def testMapConstructor() {

        when:
        def wrapper = new BashWrapperBuilder(
                input: 'alpha',
                scratch: '$var_x',
                workDir: Paths.get('a'),
                targetDir: Paths.get('b'),
                container: 'docker_x',
                environment: [a:1, b:2],
                script: 'echo ciao',
                shell: 'bash -e'
        )

        then:
        wrapper.scratch == '$var_x'
        wrapper.input == 'alpha'
        wrapper.workDir == Paths.get('a')
        wrapper.targetDir == Paths.get('b')
        wrapper.dockerImage == 'docker_x'
        wrapper.environment ==  [a:1, b:2]
        wrapper.script ==  'echo ciao'
        wrapper.shell == 'bash -e'
    }


    def testBashWrapperTest () {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def bash = new BashWrapperBuilder(workDir: folder, script: 'echo Hello world!')
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))

        folder.resolve('.command.sh').text ==
                '''
                #!/bin/bash -ue
                echo Hello world!
                '''
                .stripIndent().leftTrim()


        folder.resolve('.command.run').text ==
                """
                #!/bin/bash -Eeu
                trap on_exit 1 2 3 15 ERR TERM USR1 USR2
                function on_exit() { local exit_status=\${1:-\$?}; printf \$exit_status > ${folder}/.exitcode; exit \$exit_status; }
                touch ${folder}/.command.begin
                [ -f ${folder}/.command.env ] && source ${folder}/.command.env
                ( /bin/bash -ue ${folder}/.command.sh ) &> ${folder}/.command.out
                on_exit
                """
                .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }



    def testBashWithDockerTest () {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * bash run through docker
         */
        when:
        def bash = new BashWrapperBuilder(
                name: 'xyz',
                workDir: folder,
                script: 'echo Hello world!',
                dockerMount: Paths.get('/some/path'),
                dockerConfig: [image: 'busybox', temp: 'auto', sudo: true, enabled: true] )
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))

        folder.resolve('.command.sh').text ==
                '''
                #!/bin/bash -ue
                echo Hello world!
                '''
                        .stripIndent().leftTrim()


        folder.resolve('.command.run').text ==
                """
                #!/bin/bash -Eeu
                trap on_exit 1 2 3 15 ERR TERM USR1 USR2
                function on_exit() { local exit_status=\${1:-\$?}; printf \$exit_status > ${folder}/.exitcode; exit \$exit_status; }
                touch ${folder}/.command.begin
                ( sudo docker run -v \$(mktemp -d):/tmp -v /some/path:/some/path -v \$PWD:\$PWD -w \$PWD --name xyz busybox /bin/bash -ue ${folder}/.command.sh ) &> ${folder}/.command.out
                sudo docker rm xyz &>/dev/null || true &
                on_exit
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    def testBashWithDockerTest2() {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * bash run through docker
         */
        when:
        def bash = new BashWrapperBuilder(
                name: 'xyz',
                workDir: folder,
                script: 'echo Hello world!',
                dockerMount: Paths.get('/some/path'),
                dockerConfig: [image: 'busybox', temp: 'auto', enabled: true] )
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))

        folder.resolve('.command.sh').text ==
                '''
                #!/bin/bash -ue
                echo Hello world!
                '''
                        .stripIndent().leftTrim()


        folder.resolve('.command.run').text ==
                """
                #!/bin/bash -Eeu
                trap on_exit 1 2 3 15 ERR TERM USR1 USR2
                function on_exit() { local exit_status=\${1:-\$?}; printf \$exit_status > ${folder}/.exitcode; exit \$exit_status; }
                touch ${folder}/.command.begin
                ( docker run -v \$(mktemp -d):/tmp -v /some/path:/some/path -v \$PWD:\$PWD -w \$PWD --name xyz busybox /bin/bash -ue ${folder}/.command.sh ) &> ${folder}/.command.out
                docker rm xyz &>/dev/null || true &
                on_exit
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    def testBashWithDockerTest3() {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * bash run through docker
         */
        when:
        def bash = new BashWrapperBuilder(
                name: 'c1',
                workDir: folder,
                script: 'echo Hello world!',
                dockerMount: Paths.get('/some/path'),
                dockerConfig: [image: 'ubuntu', temp: 'auto', enabled: true, remove:false] )
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))

        folder.resolve('.command.sh').text ==
                '''
                #!/bin/bash -ue
                echo Hello world!
                '''
                        .stripIndent().leftTrim()


        folder.resolve('.command.run').text ==
                """
                #!/bin/bash -Eeu
                trap on_exit 1 2 3 15 ERR TERM USR1 USR2
                function on_exit() { local exit_status=\${1:-\$?}; printf \$exit_status > ${folder}/.exitcode; exit \$exit_status; }
                touch ${folder}/.command.begin
                ( docker run -v \$(mktemp -d):/tmp -v /some/path:/some/path -v \$PWD:\$PWD -w \$PWD --name c1 ubuntu /bin/bash -ue ${folder}/.command.sh ) &> ${folder}/.command.out
                on_exit
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

}
