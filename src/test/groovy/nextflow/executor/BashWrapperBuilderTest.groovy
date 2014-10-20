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

    def 'test change to scratchDir' () {

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


    def 'test map constructor'() {

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


    def 'test bash wrapper' () {

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
                #!/bin/bash -ue
                nxf_kill() {
                    declare -a ALL_CHILD
                    while read P PP;do
                        ALL_CHILD[\$PP]+=" \$P"
                    done < <(ps -e -o pid= -o ppid=)

                    walk() {
                        [[ \$1 != \$\$ ]] && kill \$1 2>/dev/null || true
                        for i in \${ALL_CHILD[\$1]:=}; do walk \$i; done
                    }

                    walk \$1
                }

                on_exit() {
                  set +e
                  local exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    [[ "\$pid" ]] && nxf_kill \$pid
                }

                trap on_exit EXIT
                trap on_term TERM USR1 USR2

                touch ${folder}/.command.begin
                [ -f ${folder}/.command.env ] && source ${folder}/.command.env

                set +e
                (
                /bin/bash -ue ${folder}/.command.sh &> ${folder}/.command.out
                ) &
                pid=\$!
                wait \$pid || ret=\$?
                """
                .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    def 'test bash wrapper with trace'() {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def bash = new BashWrapperBuilder(workDir: folder, script: 'echo Hello world!', statsEnabled: true)
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
                 #!/bin/bash -ue
                nxf_kill() {
                    declare -a ALL_CHILD
                    while read P PP;do
                        ALL_CHILD[\$PP]+=" \$P"
                    done < <(ps -e -o pid= -o ppid=)

                    walk() {
                        [[ \$1 != \$\$ ]] && kill \$1 2>/dev/null || true
                        for i in \${ALL_CHILD[\$1]:=}; do walk \$i; done
                    }

                    walk \$1
                }

                on_exit() {
                  set +e
                  local exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    [[ "\$pid" ]] && nxf_kill \$pid
                }

                trap on_exit EXIT
                trap on_term TERM USR1 USR2

                touch ${folder}/.command.begin
                [ -f ${folder}/.command.env ] && source ${folder}/.command.env

                set +e
                (
                /bin/bash -ue ${folder}/.command.run.1
                ) &
                pid=\$!
                wait \$pid || ret=\$?
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    /**
     * test running with Docker executed as 'sudo'
     */
    def 'test bash wrapper with docker' () {

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
                #!/bin/bash -ue
                nxf_kill() {
                    declare -a ALL_CHILD
                    while read P PP;do
                        ALL_CHILD[\$PP]+=" \$P"
                    done < <(ps -e -o pid= -o ppid=)

                    walk() {
                        [[ \$1 != \$\$ ]] && kill \$1 2>/dev/null || true
                        for i in \${ALL_CHILD[\$1]:=}; do walk \$i; done
                    }

                    walk \$1
                }

                on_exit() {
                  set +e
                  local exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    sudo docker stop -t 1 xyz
                }

                trap on_exit EXIT
                trap on_term TERM USR1 USR2

                touch ${folder}/.command.begin

                set +e
                (
                sudo docker run -v \$(mktemp -d):/tmp -v /some/path:/some/path -v \$PWD:\$PWD -w \$PWD --name xyz busybox /bin/bash -ue ${folder}/.command.sh &> ${folder}/.command.out
                ) &
                pid=\$!
                wait \$pid || ret=\$?
                sudo docker rm xyz &>/dev/null &
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    def 'test bash wrapper with docker 2'() {

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
                #!/bin/bash -ue
                nxf_kill() {
                    declare -a ALL_CHILD
                    while read P PP;do
                        ALL_CHILD[\$PP]+=" \$P"
                    done < <(ps -e -o pid= -o ppid=)

                    walk() {
                        [[ \$1 != \$\$ ]] && kill \$1 2>/dev/null || true
                        for i in \${ALL_CHILD[\$1]:=}; do walk \$i; done
                    }

                    walk \$1
                }

                on_exit() {
                  set +e
                  local exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    docker stop -t 1 xyz
                }

                trap on_exit EXIT
                trap on_term TERM USR1 USR2

                touch ${folder}/.command.begin

                set +e
                (
                docker run -v \$(mktemp -d):/tmp -v /some/path:/some/path -v \$PWD:\$PWD -w \$PWD --name xyz busybox /bin/bash -ue ${folder}/.command.sh &> ${folder}/.command.out
                ) &
                pid=\$!
                wait \$pid || ret=\$?
                docker rm xyz &>/dev/null &
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    /**
     * Test run in a docker container, without removing it
     */
    def 'test bash wrapper with docker 3'() {

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
                #!/bin/bash -ue
                nxf_kill() {
                    declare -a ALL_CHILD
                    while read P PP;do
                        ALL_CHILD[\$PP]+=" \$P"
                    done < <(ps -e -o pid= -o ppid=)

                    walk() {
                        [[ \$1 != \$\$ ]] && kill \$1 2>/dev/null || true
                        for i in \${ALL_CHILD[\$1]:=}; do walk \$i; done
                    }

                    walk \$1
                }

                on_exit() {
                  set +e
                  local exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    docker stop -t 1 c1
                }

                trap on_exit EXIT
                trap on_term TERM USR1 USR2

                touch ${folder}/.command.begin

                set +e
                (
                docker run -v \$(mktemp -d):/tmp -v /some/path:/some/path -v \$PWD:\$PWD -w \$PWD --name c1 ubuntu /bin/bash -ue ${folder}/.command.sh &> ${folder}/.command.out
                ) &
                pid=\$!
                wait \$pid || ret=\$?
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    def 'test shell exit function' () {

        def bash

        when:
        bash = [:] as BashWrapperBuilder
        then:
        bash.scriptCleanUp( Paths.get('/my/exit/file'), null ) ==
                    '''
                    nxf_kill() {
                        declare -a ALL_CHILD
                        while read P PP;do
                            ALL_CHILD[$PP]+=" $P"
                        done < <(ps -e -o pid= -o ppid=)

                        walk() {
                            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
                            for i in ${ALL_CHILD[$1]:=}; do walk $i; done
                        }

                        walk $1
                    }

                    on_exit() {
                      set +e
                      local exit_status=${ret:=$?}
                      printf $exit_status > /my/exit/file
                      exit $exit_status
                    }

                    on_term() {
                        set +e
                        [[ "$pid" ]] && nxf_kill $pid
                    }

                    trap on_exit EXIT
                    trap on_term TERM USR1 USR2
                    '''
                    .stripIndent().leftTrim()


        when:
        bash = [:] as BashWrapperBuilder
        then:
        bash.scriptCleanUp( Paths.get('/my/exit/xxx'), 'docker stop x' ) ==
                '''
                    nxf_kill() {
                        declare -a ALL_CHILD
                        while read P PP;do
                            ALL_CHILD[$PP]+=" $P"
                        done < <(ps -e -o pid= -o ppid=)

                        walk() {
                            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
                            for i in ${ALL_CHILD[$1]:=}; do walk $i; done
                        }

                        walk $1
                    }

                    on_exit() {
                      set +e
                      local exit_status=${ret:=$?}
                      printf $exit_status > /my/exit/xxx
                      exit $exit_status
                    }

                    on_term() {
                        set +e
                        docker stop x
                    }

                    trap on_exit EXIT
                    trap on_term TERM USR1 USR2
                    '''
                        .stripIndent().leftTrim()

    }

}
