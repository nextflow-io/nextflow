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
                trap on_term TERM INT USR1 USR2

                touch ${folder}/.command.begin
                [ -f ${folder}/.command.env ] && source ${folder}/.command.env

                set +e
                (
                /bin/bash -ue ${folder}/.command.sh &> .command.out
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
                trap on_term TERM INT USR1 USR2

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

    def 'test bash wrapper with scratch and input'() {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def bash = new BashWrapperBuilder(workDir: folder, script: 'echo Hello world!', scratch: true, input: 'Ciao ciao')
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))
        Files.exists(folder.resolve('.command.in'))

        folder.resolve('.command.in').text == 'Ciao ciao'

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
                trap on_term TERM INT USR1 USR2

                touch ${folder}/.command.begin
                [ -f ${folder}/.command.env ] && source ${folder}/.command.env
                NXF_SCRATCH=\${TMPDIR:-`mktemp -d`} && cd \$NXF_SCRATCH

                set +e
                (
                /bin/bash -ue ${folder}/.command.sh < ${folder}/.command.in &> .command.out
                ) &
                pid=\$!
                wait \$pid || ret=\$?
                cp .command.out ${folder} || true
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    def 'test bash wrapper with scratch and input and stats'() {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def bash = new BashWrapperBuilder(
                workDir: folder,
                script: 'echo Hello world!',
                scratch: true,
                input: 'data xyz',
                statsEnabled: true
        )
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))
        Files.exists(folder.resolve('.command.in'))
        Files.exists(folder.resolve('.command.run.1'))

        /*
         * data input file
         */
        folder.resolve('.command.in').text == 'data xyz'

        /*
         * the user scritp  file
         */
        folder.resolve('.command.sh').text ==
                '''
                #!/bin/bash -ue
                echo Hello world!
                '''
                        .stripIndent().leftTrim()


        /*
         * the main script launcher
         */
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
                trap on_term TERM INT USR1 USR2

                touch ${folder}/.command.begin
                [ -f ${folder}/.command.env ] && source ${folder}/.command.env
                NXF_SCRATCH=\${TMPDIR:-`mktemp -d`} && cd \$NXF_SCRATCH

                set +e
                (
                /bin/bash -ue ${folder}/.command.run.1
                ) &
                pid=\$!
                wait \$pid || ret=\$?
                cp .command.out ${folder} || true
                cp .command.trace ${folder} || true
                """
                        .stripIndent().leftTrim()

        folder.resolve('.command.run.1').text ==
            """
            #!/bin/bash -ue
            nxf_tree() {
                declare -a ALL_CHILD
                while read P PP;do
                    ALL_CHILD[\$PP]+=" \$P"
                done < <(ps -e -o pid= -o ppid=)

                stat() {
                    local x_ps=\$(ps -o pid=,state=,pcpu=,pmem=,vsz=,rss= \$1)
                    local x_io=\$(cat /proc/\$1/io 2> /dev/null | sed 's/^.*:\\s*//' | tr '\\n' ' ')
                    local x_vm=\$(cat /proc/\$1/status 2> /dev/null | egrep 'VmPeak|VmHWM' | sed 's/^.*:\\s*//' | sed 's/[\\sa-zA-Z]*\$//' | tr '\\n' ' ')
                    [[ ! \$x_ps ]] && return 0

                    printf "\$x_ps"
                    if [[ \$x_vm ]]; then printf " \$x_vm"; else printf " 0 0"; fi
                    if [[ \$x_io ]]; then printf " \$x_io"; fi
                    printf "\\n"
                }

                walk() {
                    stat \$1
                    for i in \${ALL_CHILD[\$1]:=}; do walk \$i; done
                }

                walk \$1
            }

            nxf_pstat() {
                local data=\$(nxf_tree \$1)
                local tot=''
                if [[ "\$data" ]]; then
                  tot=\$(awk '{ t3+=(\$3*10); t4+=(\$4*10); t5+=\$5; t6+=\$6; t7+=\$7; t8+=\$8; t9+=\$9; t10+=\$10; t11+=\$11; t12+=\$12; t13+=\$13; t14+=\$14 } END { print NR,"0",t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14 }' <<< "\$data")
                  printf "\$tot\\n"
                fi
            }

            nxf_sleep() {
              if [[ \$1 < 0 ]]; then sleep 5;
              elif [[ \$1 < 10 ]]; then sleep 0.1;
              elif [[ \$1 < 130 ]]; then sleep 1;
              else sleep 5; fi
            }

            nxf_date() {
                case `uname` in
                    Darwin) if hash gdate 2>/dev/null; then echo 'gdate +%s%3N'; else echo 'date +%s000'; fi;;
                    *) echo 'date +%s%3N';;
                esac
            }

            NXF_DATE=\$(nxf_date)

            nxf_trace() {
              local pid=\$1; local trg=\$2;
              local tot;
              local count=0;
              declare -a max=(); for i in {0..13}; do max[i]=0; done
              while [[ true ]]; do
                tot=\$(nxf_pstat \$pid)
                [[ ! \$tot ]] && break
                IFS=' ' read -a val <<< "\$tot"; unset IFS
                for i in {0..13}; do
                  [ \${val[i]} -gt \${max[i]} ] && max[i]=\${val[i]}
                done
                echo "pid state %cpu %mem vmem rss peak_vmem peak_rss rchar wchar syscr syscw read_bytes write_bytes" > \$trg
                echo "\${max[@]}" >> \$trg
                nxf_sleep \$count
                count=\$((count+1))
              done
            }


            trap 'exit \${ret:=\$?}' EXIT
            start_millis=\$(\$NXF_DATE)
            (
            /bin/bash -ue ${folder}/.command.sh < ${folder}/.command.in &> .command.out
            ) &
            pid=\$!
            nxf_trace "\$pid" .command.trace &
            mon=\$!
            wait \$pid
            ret=\$?
            end_millis=\$(\$NXF_DATE)
            kill \$mon || wait \$mon
            [ -f .command.trace ] && echo \$((end_millis-start_millis)) >> .command.trace
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
                    sudo docker kill xyz
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                touch ${folder}/.command.begin

                set +e
                (
                sudo docker run -i -v \$(mktemp -d):/tmp -v /some/path:/some/path -v \$PWD:\$PWD -w \$PWD --name xyz busybox /bin/bash -ue ${folder}/.command.sh &> .command.out
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
                    docker kill xyz
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                touch ${folder}/.command.begin

                set +e
                (
                docker run -i -v \$(mktemp -d):/tmp -v /some/path:/some/path -v \$PWD:\$PWD -w \$PWD --name xyz busybox /bin/bash -ue ${folder}/.command.sh &> .command.out
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
                    docker kill c1
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                touch ${folder}/.command.begin

                set +e
                (
                docker run -i -v \$(mktemp -d):/tmp -v /some/path:/some/path -v \$PWD:\$PWD -w \$PWD --name c1 ubuntu /bin/bash -ue ${folder}/.command.sh &> .command.out
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
                    trap on_term TERM INT USR1 USR2
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
                    trap on_term TERM INT USR1 USR2
                    '''
                        .stripIndent().leftTrim()

    }

}
