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
                #!/bin/bash
                set -e
                function on_exit() {
                  set +e
                  local exit_status
                  local err; err=\${1:-1} && [[ \$err == 0 ]] && err=1
                  local exit_file=\${2:-.exitfile}
                  if [[ \$ret ]];
                  then exit_status=\$ret;
                  else exit_status=\$err; fi
                  [[ \$pid ]] && kill -0 \$pid &>/dev/null && kill \$pid
                  printf \$exit_status > \$exit_file; exit \$exit_status;
                }

                trap 'on_exit \$? ${folder}/.exitcode' EXIT
                touch ${folder}/.command.begin
                [ -f ${folder}/.command.env ] && source ${folder}/.command.env

                # Launch job execution -- null
                ( /bin/bash -ue ${folder}/.command.sh &> ${folder}/.command.out ) &
                pid=\$!

                # Finalization
                if [[ \$pid ]]; then wait \$pid; ret=\$?; fi
                """
                .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    def testBashWrapperWithStats() {

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
                #!/bin/bash
                set -e
                function on_exit() {
                  set +e
                  local exit_status
                  local err; err=\${1:-1} && [[ \$err == 0 ]] && err=1
                  local exit_file=\${2:-.exitfile}
                  if [[ \$ret ]];
                  then exit_status=\$ret;
                  else exit_status=\$err; fi
                  [[ \$pid ]] && kill -0 \$pid &>/dev/null && kill \$pid
                  printf \$exit_status > \$exit_file; exit \$exit_status;
                }

                trap 'on_exit \$? ${folder}/.exitcode' EXIT
                touch ${folder}/.command.begin
                [ -f ${folder}/.command.env ] && source ${folder}/.command.env

                # Launch job execution -- null
                ( /bin/bash -ue ${folder}/.command.sh &> ${folder}/.command.out ) &
                pid=\$!

                # Collect proc stats
                function on_trace() {
                  local pid=\$1; local trg=\$2; local xio; local xst; local temp;
                  declare -a max=(0 0 0 0)
                  while [[ \$pid ]]; do
                    local xps; xps="\$(ps -o pcpu=,pmem=,rss=,vsz=,state= \$pid)"
                    [ \$? -ne 0 ] && break
                    IFS=' ' read -a val <<< "\$xps"; unset IFS
                    for i in {0..3}; do [ \$(echo "\${val[i]} > \${max[i]}" |bc) -eq 1 ] && max[i]=\${val[i]}; done
                    temp="\$(cat /proc/\$pid/io 2> /dev/null)" && xio=\$temp
                    temp="\$(cat /proc/\$pid/status 2> /dev/null)" && xst=\$temp
                    echo -e "%cpu %mem rss vmem state" > \$trg
                    echo -e "\${val[@]}" >> \$trg
                    echo -e "\${max[@]}" >> \$trg
                    [[ \$xio ]] && echo -e "\$xio" >> \$trg
                    [[ \$xst ]] && echo -e "\$xst" >> \$trg
                    sleep 2
                  done
                }

                ( on_trace \$pid .command.trace &> /dev/null ) &
                disown

                # Finalization
                if [[ \$pid ]]; then wait \$pid; ret=\$?; fi
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    /**
     * test running with Docker executed as 'sudo'
     */
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
                #!/bin/bash
                set -e
                function on_exit() {
                  set +e
                  local exit_status
                  local err; err=\${1:-1} && [[ \$err == 0 ]] && err=1
                  local exit_file=\${2:-.exitfile}
                  if [[ \$ret ]];
                  then exit_status=\$ret;
                  else exit_status=\$err; fi
                  [[ \$pid ]] && kill -0 \$pid &>/dev/null && kill \$pid
                  printf \$exit_status > \$exit_file; exit \$exit_status;
                }

                trap 'on_exit \$? ${folder}/.exitcode' EXIT
                touch ${folder}/.command.begin

                # Launch job execution -- xyz
                ( sudo docker run -v \$(mktemp -d):/tmp -v /some/path:/some/path -v \$PWD:\$PWD -w \$PWD --name xyz busybox /bin/bash -ue ${folder}/.command.sh &> ${folder}/.command.out ) &
                pid=\$!

                # Finalization
                if [[ \$pid ]]; then wait \$pid; ret=\$?; fi
                sudo docker rm xyz &>/dev/null || true &
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
                #!/bin/bash
                set -e
                function on_exit() {
                  set +e
                  local exit_status
                  local err; err=\${1:-1} && [[ \$err == 0 ]] && err=1
                  local exit_file=\${2:-.exitfile}
                  if [[ \$ret ]];
                  then exit_status=\$ret;
                  else exit_status=\$err; fi
                  [[ \$pid ]] && kill -0 \$pid &>/dev/null && kill \$pid
                  printf \$exit_status > \$exit_file; exit \$exit_status;
                }

                trap 'on_exit \$? ${folder}/.exitcode' EXIT
                touch ${folder}/.command.begin

                # Launch job execution -- xyz
                ( docker run -v \$(mktemp -d):/tmp -v /some/path:/some/path -v \$PWD:\$PWD -w \$PWD --name xyz busybox /bin/bash -ue ${folder}/.command.sh &> ${folder}/.command.out ) &
                pid=\$!

                # Finalization
                if [[ \$pid ]]; then wait \$pid; ret=\$?; fi
                docker rm xyz &>/dev/null || true &
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    /**
     * Test run in a docker container, without removing it
     */
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
                #!/bin/bash
                set -e
                function on_exit() {
                  set +e
                  local exit_status
                  local err; err=\${1:-1} && [[ \$err == 0 ]] && err=1
                  local exit_file=\${2:-.exitfile}
                  if [[ \$ret ]];
                  then exit_status=\$ret;
                  else exit_status=\$err; fi
                  [[ \$pid ]] && kill -0 \$pid &>/dev/null && kill \$pid
                  printf \$exit_status > \$exit_file; exit \$exit_status;
                }

                trap 'on_exit \$? ${folder}/.exitcode' EXIT
                touch ${folder}/.command.begin

                # Launch job execution -- c1
                ( docker run -v \$(mktemp -d):/tmp -v /some/path:/some/path -v \$PWD:\$PWD -w \$PWD --name c1 ubuntu /bin/bash -ue ${folder}/.command.sh &> ${folder}/.command.out ) &
                pid=\$!

                # Finalization
                if [[ \$pid ]]; then wait \$pid; ret=\$?; fi
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

}
