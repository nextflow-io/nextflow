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

package nextflow.executor

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import nextflow.Session
import nextflow.cloud.aws.batch.AwsOptions
import nextflow.container.ContainerConfig
import nextflow.container.DockerBuilder
import nextflow.container.SingularityBuilder
import nextflow.processor.TaskBean
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
import spock.lang.Unroll
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BashWrapperBuilderTest extends Specification {

    def 'test change to scratchDir' () {

        setup:
        def builder = new BashWrapperBuilder(new TaskBean())

        expect:
        builder.getScratchDirectoryCommand() == null

        when:
        builder.scratch = true
        then:
        builder.getScratchDirectoryCommand() == 'NXF_SCRATCH="$(set +u; nxf_mktemp $TMPDIR)"'

        when:
        builder.scratch = '$SOME_DIR'
        then:
        builder.getScratchDirectoryCommand() == 'NXF_SCRATCH="$(set +u; nxf_mktemp $SOME_DIR)"'

        when:
        builder.scratch = '/my/temp'
        then:
        builder.getScratchDirectoryCommand() == 'NXF_SCRATCH="$(set +u; nxf_mktemp /my/temp)"'

        when:
        builder.scratch = 'ram-disk'
        then:
        builder.getScratchDirectoryCommand() == 'NXF_SCRATCH="$(nxf_mktemp /dev/shm)"'

    }


    def 'test map constructor'() {

        given:
        def bean = new TaskBean(
                input: 'alpha',
                scratch: '$var_x',
                workDir: Paths.get('a'),
                targetDir: Paths.get('b'),
                containerImage: 'docker_x',
                environment: [a:1, b:2],
                script: 'echo ciao',
                shell: ['bash','-e']
        )

        when:
        def wrapper = new BashWrapperBuilder(bean)

        then:
        wrapper.scratch == '$var_x'
        wrapper.input == 'alpha'
        wrapper.workDir == Paths.get('a')
        wrapper.targetDir == Paths.get('b')
        wrapper.containerImage == 'docker_x'
        wrapper.environment ==  [a:1, b:2]
        wrapper.script ==  'echo ciao'
        wrapper.shell == ['bash','-e']
    }


    def 'test bash wrapper' () {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def bash = new BashWrapperBuilder([
                name: 'Hello 1',
                workDir: folder,
                script: 'echo Hello world!',
                headerScript: '#BSUB -x 1\n#BSUB -y 2'
            ] as TaskBean )
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
                #BSUB -x 1
                #BSUB -y 2
                # NEXTFLOW TASK: Hello 1
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    [[ "\$pid" ]] && nxf_kill \$pid
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                NXF_SCRATCH=''
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                /bin/bash -ue ${folder}/.command.sh
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                """
                .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    def 'test bash wrapper with inputs and outputs' () {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def bash = new BashWrapperBuilder([
                name: 'Hello 1',
                workDir: folder,
                targetDir: folder,
                scratch: true,
                inputFiles: ['sample_1.fq':Paths.get('/some/data/sample_1.fq'), 'sample_2.fq':Paths.get('/some/data/sample_2.fq'), ],
                outputFiles: ['test.bam','test.bai'],
                script: 'echo Hello world!',
                headerScript: '#BSUB -x 1\n#BSUB -y 2'
        ] as TaskBean )
        bash.build()

        then:
        folder.resolve('.command.run').text ==
                """
                #!/bin/bash
                #BSUB -x 1
                #BSUB -y 2
                # NEXTFLOW TASK: Hello 1
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  rm -rf \$NXF_SCRATCH || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    [[ "\$pid" ]] && nxf_kill \$pid
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                NXF_SCRATCH="\$(set +u; nxf_mktemp \$TMPDIR)"
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH
                # stage input files
                rm -f sample_1.fq
                rm -f sample_2.fq
                ln -s /some/data/sample_1.fq sample_1.fq
                ln -s /some/data/sample_2.fq sample_2.fq

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                /bin/bash -ue ${folder}/.command.sh
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                cp .command.out ${folder}/.command.out || true
                cp .command.err ${folder}/.command.err || true
                # copies output files to target
                if [[ \${ret:=0} == 0 ]]; then
                  mkdir -p ${folder}
                  cp -fRL test.bam ${folder} || true
                  cp -fRL test.bai ${folder} || true
                fi
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }


    def 'test bash wrapper with inputs and outputs with s3 target' () {

        given:
        def folder = Files.createTempDirectory('test')
        def target = Mock(Path)
        target.toString() >> '/some/bucket'

        def bean = new TaskBean([
                name: 'Hello 1',
                workDir: folder,
                targetDir: target,
                scratch: true,
                outputFiles: ['test.bam','test.bai'],
                script: 'echo Hello world!',
                headerScript: '#BSUB -x 1\n#BSUB -y 2'
            ])

        def copy = Spy(SimpleFileCopyStrategy, constructorArgs:[bean])
        copy.getPathScheme(target) >> 's3'
        copy.getAwsOptions() >> new AwsOptions()

        /*
         * simple bash run
         */
        when:
        def bash = new BashWrapperBuilder(bean,copy)
        bash.build()

        then:
        folder.resolve('.command.run').text ==
                """
                #!/bin/bash
                #BSUB -x 1
                #BSUB -y 2
                # NEXTFLOW TASK: Hello 1
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  rm -rf \$NXF_SCRATCH || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    [[ "\$pid" ]] && nxf_kill \$pid
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                NXF_SCRATCH="\$(set +u; nxf_mktemp \$TMPDIR)"
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                # aws helper
                nxf_s3_upload() {
                    local pattern=\$1
                    local s3path=\$2
                    IFS=\$'\\n'
                    for name in \$(eval "ls -1d \$pattern");do
                      if [[ -d "\$name" ]]; then
                        aws s3 cp --only-show-errors --recursive --storage-class STANDARD "\$name" "\$s3path/\$name"
                      else
                        aws s3 cp --only-show-errors --storage-class STANDARD "\$name" "\$s3path/\$name"
                      fi
                    done
                    unset IFS
                }

                nxf_s3_download() {
                    local source=\$1
                    local target=\$2
                    local file_name=\$(basename \$1)
                    local is_dir=\$(aws s3 ls \$source | grep -F "PRE \${file_name}/" -c)
                    if [[ \$is_dir == 1 ]]; then
                        aws s3 cp --only-show-errors --recursive "\$source" "\$target"
                    else 
                        aws s3 cp --only-show-errors "\$source" "\$target"
                    fi
                }
                
                nxf_parallel() {
                    local cmd=("\$@")
                    local cpus=\$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
                    local max=\$(if (( cpus>16 )); then echo 16; else echo \$cpus; fi)
                    local i=0
                    local pid=()
                    (
                    set +u
                    while ((i<\${#cmd[@]})); do
                        local copy=()
                        for x in "\${pid[@]}"; do
                          [[ -e /proc/\$x ]] && copy+=(\$x) 
                        done
                        pid=("\${copy[@]}")
                
                        if ((\${#pid[@]}>=\$max)); then 
                          sleep 1 
                        else 
                          eval "\${cmd[\$i]}" &
                          pid+=(\$!)
                          ((i+=1))
                        fi 
                    done
                    ((\${#pid[@]}>0)) && wait \${pid[@]}
                    )
                }     
                
                touch ${folder}/.command.begin
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                /bin/bash -ue ${folder}/.command.sh
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                cp .command.out ${folder}/.command.out || true
                cp .command.err ${folder}/.command.err || true
                # copies output files to target
                if [[ \${ret:=0} == 0 ]]; then
                  nxf_s3_upload 'test.bam' s3://some/bucket || true
                  nxf_s3_upload 'test.bai' s3://some/bucket || true
                fi
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
        def bash = new BashWrapperBuilder([
                name: 'Hello 1',
                workDir: folder,
                script: 'echo Hello world!',
                statsEnabled: true] as TaskBean)
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))
        Files.exists(folder.resolve('.command.stub'))

        folder.resolve('.command.sh').text ==
                '''
                #!/bin/bash -ue
                echo Hello world!
                '''
                        .stripIndent().leftTrim()


        folder.resolve('.command.run').text ==
                """
                #!/bin/bash
                # NEXTFLOW TASK: Hello 1
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    [[ "\$pid" ]] && nxf_kill \$pid
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                NXF_SCRATCH=''
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH
                        
                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                /bin/bash ${folder}/.command.stub
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
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
        def bash = new BashWrapperBuilder([
                name: 'Hello 1',
                workDir: folder,
                script: 'echo Hello world!',
                scratch: true,
                input: 'Ciao ciao'
            ] as TaskBean)
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
                #!/bin/bash
                # NEXTFLOW TASK: Hello 1
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  rm -rf \$NXF_SCRATCH || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    [[ "\$pid" ]] && nxf_kill \$pid
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                NXF_SCRATCH="\$(set +u; nxf_mktemp \$TMPDIR)"
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                /bin/bash -ue ${folder}/.command.sh < ${folder}/.command.in
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                cp .command.out ${folder}/.command.out || true
                cp .command.err ${folder}/.command.err || true
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
        def bash = new BashWrapperBuilder([
                name: 'Hello 2',
                workDir: folder,
                script: 'echo Hello world!',
                scratch: true,
                input: 'data xyz',
                statsEnabled: true
                ] as TaskBean)
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))
        Files.exists(folder.resolve('.command.in'))
        Files.exists(folder.resolve('.command.stub'))

        /*
         * data input file
         */
        folder.resolve('.command.in').text == 'data xyz'

        /*
         * the user script  file
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
                #!/bin/bash
                # NEXTFLOW TASK: Hello 2
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  rm -rf \$NXF_SCRATCH || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    [[ "\$pid" ]] && nxf_kill \$pid
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                NXF_SCRATCH="\$(set +u; nxf_mktemp \$TMPDIR)"
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH
                
                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                /bin/bash ${folder}/.command.stub
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                cp .command.out ${folder}/.command.out || true
                cp .command.err ${folder}/.command.err || true
                cp .command.trace ${folder}/.command.trace || true
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    /**
     * test running with Docker executed as 'sudo'
     */
    def 'test bash wrapper with conda env' () {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def bash = new BashWrapperBuilder([
                name: 'Hello 1',
                workDir: folder,
                script: 'echo Hello world!',
                condaEnv: Paths.get('/some/conda/env/foo')
        ] as TaskBean )
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))


        folder.resolve('.command.run').text ==
                """
                #!/bin/bash
                # NEXTFLOW TASK: Hello 1
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    [[ "\$pid" ]] && nxf_kill \$pid
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                NXF_SCRATCH=''
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                set +u
                # conda environment
                source activate /some/conda/env/foo
                set -u
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                /bin/bash -ue ${folder}/.command.sh
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
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
        def bash = new BashWrapperBuilder([
                name: 'Hello 3',
                workDir: folder,
                script: 'echo Hello world!',
                containerImage: 'busybox',
                containerConfig: [engine: 'docker', sudo: true, enabled: true],
                containerEnabled: true,
                ] as TaskBean)
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
                # NEXTFLOW TASK: Hello 3
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  sudo docker rm \$NXF_BOXID &>/dev/null || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    sudo docker kill \$NXF_BOXID
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                export NXF_BOXID="nxf-\$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A)"
                NXF_SCRATCH=''
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                sudo docker run -i -v ${folder}:${folder} -v "\$PWD":"\$PWD" -w "\$PWD" --entrypoint /bin/bash --name \$NXF_BOXID busybox -c "/bin/bash -ue ${folder}/.command.sh"
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
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
        def bash = new BashWrapperBuilder([
                name: 'Hello 4',
                workDir: folder,
                script: 'echo Hello world!',
                containerImage: 'busybox',
                containerEnabled: true,
                containerConfig: [engine: 'docker', temp: 'auto', enabled: true]
                ] as TaskBean)
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
                # NEXTFLOW TASK: Hello 4
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  docker rm \$NXF_BOXID &>/dev/null || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    docker kill \$NXF_BOXID
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                export NXF_BOXID="nxf-\$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A)"
                NXF_SCRATCH=''
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                docker run -i -v \$(nxf_mktemp):/tmp -v ${folder}:${folder} -v "\$PWD":"\$PWD" -w "\$PWD" --entrypoint /bin/bash --name \$NXF_BOXID busybox -c "/bin/bash -ue ${folder}/.command.sh"
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
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
        def bash = new BashWrapperBuilder([
                name: 'Hello 5',
                workDir: folder,
                script: 'echo Hello world!',
                containerEnabled: true,
                containerImage: 'ubuntu',
                containerConfig: [engine: 'docker', temp: 'auto', enabled: true, remove:false, kill: false]
                ] as TaskBean)
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
                # NEXTFLOW TASK: Hello 5
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    [[ "\$pid" ]] && nxf_kill \$pid
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                export NXF_BOXID="nxf-\$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A)"
                NXF_SCRATCH=''
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                docker run -i -v \$(nxf_mktemp):/tmp -v ${folder}:${folder} -v "\$PWD":"\$PWD" -w "\$PWD" --entrypoint /bin/bash --name \$NXF_BOXID ubuntu -c "/bin/bash -ue ${folder}/.command.sh"
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    def 'test bash wrapper with docker 4'() {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * bash run through docker
         */
        when:
        def bash = new BashWrapperBuilder([
                name: 'Hello 6',
                workDir: folder,
                script: 'echo Hello world!',
                containerEnabled: true,
                containerImage: 'ubuntu',
                containerConfig: [engine: 'docker', temp: 'auto', enabled: true, remove:false, kill: 'SIGXXX']
                ] as TaskBean)
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
                # NEXTFLOW TASK: Hello 6
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    docker kill -s SIGXXX \$NXF_BOXID
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                export NXF_BOXID="nxf-\$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A)"
                NXF_SCRATCH=''
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                docker run -i -v \$(nxf_mktemp):/tmp -v ${folder}:${folder} -v "\$PWD":"\$PWD" -w "\$PWD" --entrypoint /bin/bash --name \$NXF_BOXID ubuntu -c "/bin/bash -ue ${folder}/.command.sh"
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    /**
     * test running with Docker executed as 'sudo'
     */
    def 'test bash wrapper with docker mount' () {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * bash run through docker
         */
        when:
        def bash = new BashWrapperBuilder([
                name: 'Hello 7',
                workDir: folder,
                script: 'echo Hello world!',
                containerEnabled: true,
                containerImage: 'busybox',
                containerMount: '/folder with blanks' as Path,
                containerConfig: [engine: 'docker', enabled: true]
        ] as TaskBean)
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
                # NEXTFLOW TASK: Hello 7
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  docker rm \$NXF_BOXID &>/dev/null || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    docker kill \$NXF_BOXID
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                export NXF_BOXID="nxf-\$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A)"
                NXF_SCRATCH=''
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                docker run -i -v /folder\\ with\\ blanks:/folder\\ with\\ blanks -v ${folder}:${folder} -v "\$PWD":"\$PWD" -w "\$PWD" --entrypoint /bin/bash --name \$NXF_BOXID busybox -c "/bin/bash -ue ${folder}/.command.sh"
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    def 'test bash wrapper with docker and scratch dir' () {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * bash run through docker
         */
        when:
        def bash = new BashWrapperBuilder([
                name: 'Hello 3',
                workDir: folder,
                scratch: true,
                script: 'echo Hello world!',
                containerEnabled: true,
                containerImage: 'busybox',
                containerConfig: [engine: 'docker', sudo: true, enabled: true]
        ] as TaskBean)
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
                # NEXTFLOW TASK: Hello 3
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  (sudo -n true && sudo rm -rf "\$NXF_SCRATCH" || rm -rf "\$NXF_SCRATCH")&>/dev/null || true
                  sudo docker rm \$NXF_BOXID &>/dev/null || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    sudo docker kill \$NXF_BOXID
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                export NXF_BOXID="nxf-\$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A)"
                NXF_SCRATCH="\$(set +u; nxf_mktemp \$TMPDIR)"
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                sudo docker run -i -v ${folder}:${folder} -v "\$PWD":"\$PWD" -w "\$PWD" --entrypoint /bin/bash --name \$NXF_BOXID busybox -c "/bin/bash -ue ${folder}/.command.sh"
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                cp .command.out ${folder}/.command.out || true
                cp .command.err ${folder}/.command.err || true
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }


    def 'should create script for docker executable container' () {
        given:
        def folder = TestHelper.createInMemTempDir()
        def session = new Session(); session.workDir = folder
        def task = new TaskRun(
                name: 'Hello',
                script: 'FOO=bar\ndocker-io/busybox --fox --baz',
                config: [container: true],
                workDir: folder )
        task.processor = Mock(TaskProcessor)
        task.processor.getProcessEnvironment() >> [:]
        task.processor.getSession() >> session
        task.processor.getConfig() >> [:]

        when:
        def bash = new BashWrapperBuilder(task)
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))

        folder.resolve('.command.sh').text ==
                """
                #!/bin/bash -ue
                FOO=bar
                docker run -i -e "FOO=bar" -v $folder:$folder -v "\$PWD":"\$PWD" -w "\$PWD" --name \$NXF_BOXID docker-io/busybox --fox --baz
                """
                .stripIndent().leftTrim()

        folder.resolve('.command.run').text ==
                """
                #!/bin/bash
                # NEXTFLOW TASK: Hello
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  docker rm \$NXF_BOXID &>/dev/null || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    docker kill \$NXF_BOXID
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                export NXF_BOXID="nxf-\$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A)"
                NXF_SCRATCH=''
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                # task environment
                nxf_taskenv() {
                cat << EOF
                export FOO="bar"
                EOF
                }

                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                /bin/bash -ue ${folder}/.command.sh
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                """
                        .stripIndent().leftTrim()
    }

    def 'should create docker run with registry' () {
        given:
        def folder = TestHelper.createInMemTempDir()
        def task = new TaskRun(
                name: 'Hello',
                script: 'my.registry.com/docker-io/busybox --fox --baz',
                config: [container: true],
                workDir: folder )
        task.processor = Mock(TaskProcessor)
        task.processor.getProcessEnvironment() >> [:]
        task.processor.getSession() >> new Session(docker: [registry: 'registry.com'])
        task.processor.getSession().workDir = folder
        task.processor.getConfig() >> [:]

        when:
        def bash = new BashWrapperBuilder(task)
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))

        folder.resolve('.command.sh').text ==
                """
                #!/bin/bash -ue
                docker run -i -v $folder:$folder -v "\$PWD":"\$PWD" -w "\$PWD" --name \$NXF_BOXID my.registry.com/docker-io/busybox --fox --baz
                """
                        .stripIndent().leftTrim()
    }

    def 'test bash wrapper with docker container custom options' () {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * bash run through docker
         */
        when:
        def bash = new BashWrapperBuilder([
                name: 'Hello 3',
                workDir: folder,
                scratch: true,
                script: 'echo Hello world!',
                containerEnabled: true,
                containerImage: 'busybox',
                containerOptions: '-v /foo:/bar',
                containerConfig: [engine: 'docker', enabled: true]
        ] as TaskBean)
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))


        folder.resolve('.command.run').text ==
                """
                #!/bin/bash
                # NEXTFLOW TASK: Hello 3
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  (sudo -n true && sudo rm -rf "\$NXF_SCRATCH" || rm -rf "\$NXF_SCRATCH")&>/dev/null || true
                  docker rm \$NXF_BOXID &>/dev/null || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    docker kill \$NXF_BOXID
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                export NXF_BOXID="nxf-\$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A)"
                NXF_SCRATCH="\$(set +u; nxf_mktemp \$TMPDIR)"
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                docker run -i -v ${folder}:${folder} -v "\$PWD":"\$PWD" -w "\$PWD" --entrypoint /bin/bash -v /foo:/bar --name \$NXF_BOXID busybox -c "/bin/bash -ue ${folder}/.command.sh"
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                cp .command.out ${folder}/.command.out || true
                cp .command.err ${folder}/.command.err || true
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    def 'test shell exit function' () {

        def bash

        when:
        bash = new BashWrapperBuilder( new TaskBean() )
        then:
        bash.scriptCleanUp( Paths.get("/my/exit/file's"), 'NXF_SCRATCH=xx' ) ==
                    '''
                    nxf_env() {
                        echo '============= task environment ============='
                        env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                        echo '============= task output =================='
                    }

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

                    nxf_mktemp() {
                        local base=\${1:-/tmp}
                        if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                        else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                        fi
                    }

                    on_exit() {
                      exit_status=${ret:=$?}
                      printf $exit_status > /my/exit/file\\'s
                      set +u
                      [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                      [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                      [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                      rm -rf $NXF_SCRATCH || true
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
        def builder = Mock(DockerBuilder)
        builder.getRemoveCommand() >> 'docker rm x'
        builder.getKillCommand() >> 'docker kill x'
        bash = new BashWrapperBuilder(new TaskBean())
        bash.containerBuilder = builder
        then:
        bash.scriptCleanUp( Paths.get('/my/exit/xxx'), null ) ==
                '''
                    nxf_env() {
                        echo '============= task environment ============='
                        env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                        echo '============= task output =================='
                    }

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

                    nxf_mktemp() {
                        local base=\${1:-/tmp}
                        if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                        else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                        fi
                    }

                    on_exit() {
                      exit_status=${ret:=$?}
                      printf $exit_status > /my/exit/xxx
                      set +u
                      [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                      [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                      [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                      docker rm x &>/dev/null || true
                      exit $exit_status
                    }

                    on_term() {
                        set +e
                        docker kill x
                    }

                    trap on_exit EXIT
                    trap on_term TERM INT USR1 USR2
                    '''
                        .stripIndent().leftTrim()


        when:
        builder = Mock(DockerBuilder)
        builder.getRemoveCommand() >> 'docker rm x'
        builder.getKillCommand() >> 'docker kill x'
        bash = new BashWrapperBuilder(new TaskBean())
        bash.containerBuilder = builder
        then:
        bash.scriptCleanUp( Paths.get('/my/exit/xxx'), 'NXF_SCRATCH=xxx' ) ==
                '''
                    nxf_env() {
                        echo '============= task environment ============='
                        env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                        echo '============= task output =================='
                    }

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

                    nxf_mktemp() {
                        local base=\${1:-/tmp}
                        if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                        else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                        fi
                    }

                    on_exit() {
                      exit_status=${ret:=$?}
                      printf $exit_status > /my/exit/xxx
                      set +u
                      [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                      [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                      [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                      (sudo -n true && sudo rm -rf "$NXF_SCRATCH" || rm -rf "$NXF_SCRATCH")&>/dev/null || true
                      docker rm x &>/dev/null || true
                      exit $exit_status
                    }

                    on_term() {
                        set +e
                        docker kill x
                    }

                    trap on_exit EXIT
                    trap on_term TERM INT USR1 USR2
                    '''
                        .stripIndent().leftTrim()

    }

    def 'test environment and modules' () {

        given:
        def folder = TestHelper.createInMemTempDir()

        when:
        new BashWrapperBuilder([
                workDir: folder,
                environment: [DELTA:1, OMEGA:2, BRAVO: 'hola'],
                script: 'Hello world',
                moduleNames: ['xx/1.2','yy/3.4'] ] as TaskBean ) .build()

        then:
        folder.resolve('.command.run').text == """
                    #!/bin/bash
                    # NEXTFLOW TASK: null
                    set -e
                    set -u
                    NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                    nxf_env() {
                        echo '============= task environment =============\'
                        env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                        echo '============= task output ==================\'
                    }

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

                    nxf_mktemp() {
                        local base=\${1:-/tmp}
                        if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                        else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                        fi
                    }

                    on_exit() {
                      exit_status=\${ret:=\$?}
                      printf \$exit_status > ${folder}/.exitcode
                      set +u
                      [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                      [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                      [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                      exit \$exit_status
                    }

                    on_term() {
                        set +e
                        [[ "\$pid" ]] && nxf_kill \$pid
                    }

                    trap on_exit EXIT
                    trap on_term TERM INT USR1 USR2

                    NXF_SCRATCH=\'\'
                    [[ \$NXF_DEBUG > 0 ]] && nxf_env
                    touch ${folder}/.command.begin
                    set +u
                    nxf_module_load(){
                      local mod=\$1
                      local ver=\${2:-}
                      local new_module="\$mod"; [[ \$ver ]] && new_module+="/\$ver"

                      if [[ ! \$(module list 2>&1 | grep -o "\$new_module") ]]; then
                        old_module=\$(module list 2>&1 | grep -Eow "\$mod\\/[^\\( \\n]+" || true)
                        if [[ \$ver && \$old_module ]]; then
                          module switch \$old_module \$new_module
                        else
                          module load \$new_module
                        fi
                      fi
                    }

                    nxf_module_load xx 1.2
                    nxf_module_load yy 3.4

                    set -u
                    # task environment
                    export DELTA="1"
                    export OMEGA="2"
                    export BRAVO="hola"

                    [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                    set +e
                    ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                    cout=\$ctmp/.command.out; mkfifo \$cout
                    cerr=\$ctmp/.command.err; mkfifo \$cerr
                    tee .command.out < \$cout &
                    tee1=\$!
                    tee .command.err < \$cerr >&2 &
                    tee2=\$!
                    (
                    /bin/bash -ue ${folder}/.command.sh
                    ) >\$cout 2>\$cerr &
                    pid=\$!
                    wait \$pid || ret=\$?
                    wait \$tee1 \$tee2
                    """
                .stripIndent().leftTrim()

        when:
        folder = TestHelper.createInMemTempDir()
        new BashWrapperBuilder([
                workDir: folder,
                script: 'Hello world',
                moduleNames: ['ciao/1','mondo/2', 'bioinfo-tools']
        ] as TaskBean) .build()

        then:
        folder.resolve('.command.run').text == """
                    #!/bin/bash
                    # NEXTFLOW TASK: null
                    set -e
                    set -u
                    NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                    nxf_env() {
                        echo '============= task environment =============\'
                        env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                        echo '============= task output ==================\'
                    }

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

                    nxf_mktemp() {
                        local base=\${1:-/tmp}
                        if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                        else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                        fi
                    }

                    on_exit() {
                      exit_status=\${ret:=\$?}
                      printf \$exit_status > ${folder}/.exitcode
                      set +u
                      [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                      [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                      [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                      exit \$exit_status
                    }

                    on_term() {
                        set +e
                        [[ "\$pid" ]] && nxf_kill \$pid
                    }

                    trap on_exit EXIT
                    trap on_term TERM INT USR1 USR2

                    NXF_SCRATCH=\'\'
                    [[ \$NXF_DEBUG > 0 ]] && nxf_env
                    touch ${folder}/.command.begin
                    set +u
                    nxf_module_load(){
                      local mod=\$1
                      local ver=\${2:-}
                      local new_module="\$mod"; [[ \$ver ]] && new_module+="/\$ver"

                      if [[ ! \$(module list 2>&1 | grep -o "\$new_module") ]]; then
                        old_module=\$(module list 2>&1 | grep -Eow "\$mod\\/[^\\( \\n]+" || true)
                        if [[ \$ver && \$old_module ]]; then
                          module switch \$old_module \$new_module
                        else
                          module load \$new_module
                        fi
                      fi
                    }

                    nxf_module_load ciao 1
                    nxf_module_load mondo 2
                    nxf_module_load bioinfo-tools

                    set -u
                    [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                    set +e
                    ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                    cout=\$ctmp/.command.out; mkfifo \$cout
                    cerr=\$ctmp/.command.err; mkfifo \$cerr
                    tee .command.out < \$cout &
                    tee1=\$!
                    tee .command.err < \$cerr >&2 &
                    tee2=\$!
                    (
                    /bin/bash -ue ${folder}/.command.sh
                    ) >\$cout 2>\$cerr &
                    pid=\$!
                    wait \$pid || ret=\$?
                    wait \$tee1 \$tee2
                    """
                .stripIndent().leftTrim()

    }


    def 'test before/after script' () {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * bash run through docker
         */
        when:
        def bash = new BashWrapperBuilder([
                name: 'Hello 9',
                workDir: folder,
                script: 'echo Hello world!',
                beforeScript: "init this",
                afterScript: "cleanup that"
                ] as TaskBean)
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
                # NEXTFLOW TASK: Hello 9
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    [[ "\$pid" ]] && nxf_kill \$pid
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                NXF_SCRATCH=''
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                set +u
                # user `beforeScript`
                init this
                set -u
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                /bin/bash -ue ${folder}/.command.sh
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                # user `afterScript`
                cleanup that
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()


    }

    def 'test bash wrapper with shifter'() {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * bash run through docker
         */
        when:
        def bash = new BashWrapperBuilder([
                name: 'Hello 1',
                workDir: folder,
                script: 'echo Hello world!',
                containerEnabled: true,
                containerImage: 'docker:ubuntu:latest',
                environment: [PATH: '/path/to/bin:$PATH', FOO: 'xxx'],
                containerConfig: [enabled: true, engine: 'shifter'] as ContainerConfig
        ] as TaskBean)
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
                # NEXTFLOW TASK: Hello 1
                set -e
                set -u
                NXF_DEBUG=\${NXF_DEBUG:=0}; [[ \$NXF_DEBUG > 1 ]] && set -x


                function shifter_img() {
                  local cmd=\$1
                  local image=\$2
                  shifterimg -v \$cmd \$image |  awk -F: '\$0~/"status":/{gsub("[\\", ]","",\$2);print \$2}'
                }

                function shifter_pull() {
                  local image=\$1
                  local STATUS=\$(shifter_img lookup \$image)
                  if [[ \$STATUS != READY && \$STATUS != '' ]]; then
                    STATUS=\$(shifter_img pull \$image)
                    while [[ \$STATUS != READY && \$STATUS != FAILURE && \$STATUS != '' ]]; do
                      sleep 5
                      STATUS=\$(shifter_img pull \$image)
                    done
                  fi

                  [[ \$STATUS == FAILURE || \$STATUS == '' ]] && echo "Shifter failed to pull image \\`\$image\\`" >&2  && exit 1
                }

                nxf_env() {
                    echo '============= task environment ============='
                    env | sort | sed "s/\\(.*\\)AWS\\(.*\\)=\\(.\\{6\\}\\).*/\\1AWS\\2=\\3xxxxxxxxxxxxx/"
                    echo '============= task output =================='
                }

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

                nxf_mktemp() {
                    local base=\${1:-/tmp}
                    if [[ \$(uname) = Darwin ]]; then mktemp -d \$base/nxf.XXXXXXXXXX
                    else TMPDIR="\$base" mktemp -d -t nxf.XXXXXXXXXX
                    fi
                }

                on_exit() {
                  exit_status=\${ret:=\$?}
                  printf \$exit_status > ${folder}/.exitcode
                  set +u
                  [[ "\$tee1" ]] && kill \$tee1 2>/dev/null
                  [[ "\$tee2" ]] && kill \$tee2 2>/dev/null
                  [[ "\$ctmp" ]] && rm -rf \$ctmp || true
                  exit \$exit_status
                }

                on_term() {
                    set +e
                    [[ "\$pid" ]] && nxf_kill \$pid
                }

                trap on_exit EXIT
                trap on_term TERM INT USR1 USR2

                export NXF_BOXID="nxf-\$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A)"
                NXF_SCRATCH=''
                [[ \$NXF_DEBUG > 0 ]] && nxf_env
                touch ${folder}/.command.begin
                # task environment
                nxf_taskenv() {
                cat << EOF
                export PATH="/path/to/bin:\\\$PATH"
                export FOO="xxx"
                EOF
                }

                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH

                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                shifter_pull docker:ubuntu:latest
                shifter --image docker:ubuntu:latest /bin/bash -c "eval \$(nxf_taskenv); /bin/bash -ue ${folder}/.command.sh"
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }

    def 'should resolve foreign files' () {

        given:
        def INPUTS = [
                'foo.txt': Paths.get('/some/foo.txt'),
                'bar.txt': Paths.get('/some/bar.txt'),
        ]

        def RESOLVED = [
                'foo.txt': Paths.get('/some/foo.txt'),
                'BAR.txt': Paths.get('/some/BAR.txt'),
        ]
        def bean = Mock(TaskBean)
        bean.getInputFiles() >> INPUTS
        def copy = Mock(SimpleFileCopyStrategy)
        def bash = Spy(BashWrapperBuilder, constructorArgs:[bean, copy])

        when:
        def files = bash.getResolvedInputs()
        then:
        1 * copy.resolveForeignFiles(INPUTS) >> RESOLVED
        files == RESOLVED

    }

    def 'should create module command' () {
        given:
        def wrapper = new BashWrapperBuilder(Mock(TaskBean))

        expect:
        wrapper.moduleLoad('foo')  == 'nxf_module_load foo'
        wrapper.moduleLoad('foo/1.2')  == 'nxf_module_load foo 1.2'
        wrapper.moduleLoad('foo/bar/1.2')  == 'nxf_module_load foo/bar 1.2'
        wrapper.moduleLoad('foo/bar/')  == 'nxf_module_load foo/bar '

    }


    def 'should create the task environment' () {

        given:
        def ENV = [FOO: 'hello', BAR: 'hello world', PATH: '/some/path:$PATH']
        BashWrapperBuilder builder
        def strategy = new SimpleFileCopyStrategy()
        def env

        when:
        builder = new BashWrapperBuilder(runWithContainer: false, copyStrategy: strategy)
        env = builder.createTaskEnvironment([:])
        then:
        env == null

        when:
        builder = new BashWrapperBuilder(runWithContainer: false, copyStrategy: strategy)
        env = builder.createTaskEnvironment(ENV)
        then:
        env ==  '''
                # task environment
                export FOO="hello"
                export BAR="hello world"
                export PATH="/some/path:$PATH"
                '''
                .stripIndent().leftTrim()

        when:
        builder = new BashWrapperBuilder(runWithContainer: true, copyStrategy: strategy)
        env = builder.createTaskEnvironment([:])
        then:
        env == null

        when:
        builder = new BashWrapperBuilder(runWithContainer: true, copyStrategy: strategy)
        env = builder.createTaskEnvironment(ENV)
        then:
        env ==  '''
                # task environment
                nxf_taskenv() {
                cat << EOF
                export FOO="hello"
                export BAR="hello world"
                export PATH="/some/path:\\$PATH"
                EOF
                }
                '''
                .stripIndent().leftTrim()

    }

    @Unroll
    def 'check legacy stub script flag' () {
        given:
        def builder = Spy(BashWrapperBuilder)

        when:
        def flag = builder.isLegacyStubScript()
        then:
        builder.isMacOS() >> IS_MAC
        builder.isContainerEnabled() >> IS_CONTAINER
        flag == IS_LEGACY

        where:
        IS_MAC  | IS_CONTAINER  | IS_LEGACY
        false   | false         | false
        false   | true          | false
        true    | false         | true
        true    | true          | false
    }

    def 'should return stub script for mac' () {
        given:
        def stub = Files.createTempFile('test',null)
        def copy = Mock(ScriptFileCopyStrategy)
        def bean = Mock(TaskBean)

        def builder = Spy(BashWrapperBuilder)
        builder.copyStrategy = copy
        builder.stubFile = stub
        builder.scriptFile = Paths.get('my-command.sh')
        builder.bean = bean

        when:
        builder.createStubScript('/bin/bash')

        then:
        builder
        1 * builder.isLegacyStubScript() >> true
        1 * builder.fixOwnership() >> false
        _ * copy.fileStr(_ as Path) >>  { Path path -> path.getFileName().toString() }
        _ * copy.pipeInputFile(_ as Path) >> { " < ${it.getFileName()}"  }

        stub.text == '''
                #!/bin/bash
                set -e
                set -u
                NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 2 ]] && set -x
                
                nxf_tree() {
                    declare -a ALL_CHILD
                    while read P PP;do
                        ALL_CHILD[$PP]+=" $P"
                    done < <(ps -e -o pid= -o ppid=)
                
                    stat() {
                        local x_ps=$(ps -o pid= -o state= -o pcpu= -o pmem= -o vsz= -o rss= $1)
                        local x_io=$(cat /proc/$1/io 2> /dev/null | sed 's/^.*:\\s*//' | tr '\\n' ' ')
                        local x_vm=$(cat /proc/$1/status 2> /dev/null | egrep 'VmPeak|VmHWM' | sed 's/^.*:\\s*//' | sed 's/[\\sa-zA-Z]*$//' | tr '\\n' ' ')
                        [[ ! $x_ps ]] && return 0
                
                        printf "$x_ps"
                        if [[ $x_vm ]]; then printf " $x_vm"; else printf " 0 0"; fi
                        if [[ $x_io ]]; then printf " $x_io"; fi
                        printf "\\n"
                    }
                
                    walk() {
                        stat $1
                        for i in ${ALL_CHILD[$1]:=}; do walk $i; done
                    }
                
                    walk $1
                }
                
                nxf_pstat() {
                    local data=$(nxf_tree $1)
                    local tot=\'\'
                    if [[ "$data" ]]; then
                      tot=$(awk '{ t3+=($3*10); t4+=($4*10); t5+=$5; t6+=$6; t7+=$7; t8+=$8; t9+=$9; t10+=$10; t11+=$11; t12+=$12; t13+=$13; t14+=$14 } END { printf "%d 0 %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f\\n", NR,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14 }' <<< "$data")
                      printf "$tot\\n" || true
                    fi
                }
                
                nxf_sleep() {
                  if [[ $1 < 0 ]]; then sleep 5;
                  elif [[ $1 < 10 ]]; then sleep 0.1 2>/dev/null || sleep 1;
                  elif [[ $1 < 130 ]]; then sleep 1;
                  else sleep 5; fi
                }
                
                nxf_date() {
                    local ts=$(date +%s%3N); [[ $ts == *3N ]] && date +%s000 || echo $ts
                }
                
                nxf_trace() {
                  local pid=$1; local trg=$2;
                  local tot;
                  local count=0;
                  declare -a max=(); for i in {0..13}; do max[i]=0; done
                  while [[ true ]]; do
                    if ! kill -0 $pid 2>/dev/null; then exit 0; fi
                    tot=$(nxf_pstat $pid)
                    [[ ! $tot ]] && break
                    IFS=' ' read -a val <<< "$tot"; unset IFS
                    for i in {0..13}; do
                      [ ${val[i]} -gt ${max[i]} ] && max[i]=${val[i]}
                    done
                    echo "pid state %cpu %mem vmem rss peak_vmem peak_rss rchar wchar syscr syscw read_bytes write_bytes" > $trg
                    echo "${max[@]}" >> $trg
                    nxf_sleep $count
                    count=$((count+1))
                  done
                }
                
                
                trap 'exit ${ret:=$?}' EXIT
                touch .command.trace
                start_millis=$(nxf_date)
                (
                /bin/bash my-command.sh
                ) &
                pid=$!
                nxf_trace "$pid" .command.trace &
                mon=$!
                wait $pid || ret=$?
                end_millis=$(nxf_date)
                wait $mon
                echo $((end_millis-start_millis)) >> .command.trace
                '''
                .stripIndent().leftTrim()

        cleanup:
        stub?.delete()
    }

    def 'should return stub script for linux' () {
        given:
        def stub = Files.createTempFile('test',null)
        def copy = Mock(ScriptFileCopyStrategy)
        def bean = Mock(TaskBean)

        def builder = Spy(BashWrapperBuilder)
        builder.copyStrategy = copy
        builder.stubFile = stub
        builder.scriptFile = Paths.get('my-command.sh')
        builder.bean = bean

        when:
        builder.createStubScript('/bin/bash')

        then:
        builder
        1 * builder.isLegacyStubScript() >> false
        1 * builder.fixOwnership() >> false
        _ * copy.fileStr(_ as Path) >>  { Path it -> it.getFileName().toString() }
        _ * copy.pipeInputFile(_ as Path) >> { " < ${it.getFileName()}"  }

        stub.text == '''
                #!/bin/bash
                set -e
                set -u
                NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 2 ]] && set -x
                
                set -o pipefail
                prev_total=0
                declare -a prev_time
                mem_tot=$(< /proc/meminfo grep MemTotal | awk '{print $2}')
                num_cpus=$(< /proc/cpuinfo grep '^processor' -c)
                
                nxf_pcpu() {
                    local pid=$1
                    local proc_time=$(2> /dev/null < /proc/$pid/stat awk '{sum=$14+$15; printf "%.0f",sum}' || echo -n 0)
                    local cpu_usage=$(echo -n $proc_time ${prev_time[pid]:-0} $total_time $prev_total $num_cpus | awk '{ pct=($1-$2)/($3-$4)*$5 *100; printf "%.1f", pct }' )
                    prev_time[pid]=$proc_time
                    nxf_pcpu_ret=$cpu_usage
                }
                
                nxf_tree() {
                    declare -a ALL_CHILD
                    while read P PP;do
                        ALL_CHILD[$PP]+=" $P"
                    done < <(ps -e -o pid= -o ppid=)
                
                    stat() {
                        nxf_pcpu $1 
                        local x_pid=$1
                        local x_stat=$(2> /dev/null < /proc/$1/stat awk '{print $3}' || echo -n X)
                        local x_pcpu=$nxf_pcpu_ret
                        
                        local x_vsz=$(2> /dev/null < /proc/$1/stat awk '{printf "%.0f", $23/1024}' || echo -n 0)
                        local x_rss=$(2> /dev/null < /proc/$1/status grep VmRSS | awk '{print $2}' || echo -n 0)
                        local x_pmem=$(echo $x_rss | awk -v mem_tot=$mem_tot '{printf "%.1f", $1/mem_tot*100}')

                        local x_io=$(2> /dev/null < /proc/$1/io sed 's/^.*:\\s*//' | tr '\\n' ' ' || echo -n 0)
                        local x_vm=$(2> /dev/null < /proc/$1/status egrep 'VmPeak|VmHWM' | sed 's/^.*:\\s*//' | sed 's/[\\sa-zA-Z]*$//' | tr '\\n' ' ' || echo -n 0)
                
                        stat_ret+="$x_pid $x_stat $x_pcpu $x_pmem $x_vsz $x_rss"
                        if [[ $x_vm ]]; then stat_ret+=" $x_vm"; else stat_ret+=" 0 0"; fi
                        if [[ $x_io ]]; then stat_ret+=" $x_io"; fi
                        stat_ret+='\\n\'
                    }
                
                    walk() {
                        stat $1 
                        for i in ${ALL_CHILD[$1]:=}; do walk $i; done
                    }
                
                    stat_ret=\'\'
                    total_time=$(grep '^cpu ' /proc/stat |awk '{sum=$2+$3+$4+$5+$6+$7+$8+$9+$10; printf "%.0f",sum}')
                    walk $1
                    prev_total=$total_time
                    nxf_tree_ret=$stat_ret  
                }
                
                nxf_pstat() {
                    nxf_tree $1
                    if [[ "$nxf_tree_ret" ]]; then
                      nxf_pstat_ret=$(printf "$nxf_tree_ret" | awk '{ t3+=($3*10); t4+=($4*10); t5+=$5; t6+=$6; t7+=$7; t8+=$8; t9+=$9; t10+=$10; t11+=$11; t12+=$12; t13+=$13; t14+=$14 } END { printf "%d 0 %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f\\n", NR,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14 }')
                    else
                      nxf_pstat_ret=''  
                    fi
                }
                
                nxf_sleep() {
                  if [[ $1 < 0 ]]; then sleep 5;
                  elif [[ $1 < 10 ]]; then sleep 0.1 2>/dev/null || sleep 1;
                  elif [[ $1 < 130 ]]; then sleep 1;
                  else sleep 5; fi
                }
                
                nxf_date() {
                    local ts=$(date +%s%3N); [[ $ts == *3N ]] && date +%s000 || echo $ts
                }
                
                nxf_trace() {
                  local pid=$1; local trg=$2;
                  local count=0;
                  declare -a max=(); for i in {0..13}; do max[i]=0; done
                  while [[ -d /proc/$pid ]]; do
                    nxf_pstat $pid
                    if [[ "$nxf_pstat_ret" ]]; then
                    IFS=' ' read -a val <<< "$nxf_pstat_ret"; unset IFS
                    for i in {0..13}; do
                      [ ${val[i]} -gt ${max[i]} ] && max[i]=${val[i]}
                    done
                    echo "pid state %cpu %mem vmem rss peak_vmem peak_rss rchar wchar syscr syscw read_bytes write_bytes" > $trg
                    echo "${max[@]}" >> $trg
                    fi
                    nxf_sleep $count
                    count=$((count+1))
                  done
                } 
                
                
                trap 'exit ${ret:=$?}' EXIT
                touch .command.trace
                start_millis=$(nxf_date)
                (
                /bin/bash my-command.sh
                ) &
                pid=$!
                nxf_trace "$pid" .command.trace &
                mon=$!
                wait $pid || ret=$?
                end_millis=$(nxf_date)
                wait $mon
                echo $((end_millis-start_millis)) >> .command.trace
                '''
                .stripIndent().leftTrim()

        cleanup:
        stub?.delete()
    }

    def 'should return stub script for linux with input and fixOwnership' () {
        given:
        def stub = Files.createTempFile('test',null)
        def copy = Mock(ScriptFileCopyStrategy)
        def bean = Mock(TaskBean)

        def builder = Spy(BashWrapperBuilder)
        builder.copyStrategy = copy
        builder.stubFile = stub
        builder.scriptFile = Paths.get('my-command.sh')
        builder.bean = bean
        builder.bean.input = 'hello world'

        when:
        builder.createStubScript('/bin/bash')

        then:
        builder
        1 * builder.isLegacyStubScript() >> false
        1 * builder.fixOwnership() >> true
        _ * copy.fileStr(_ as Path) >>  { Path it -> it.getFileName().toString() }
        _ * copy.pipeInputFile(_ as Path) >> { " < ${it.getFileName()}"  }

        then:
        stub.text == '''
                #!/bin/bash
                set -e
                set -u
                NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 2 ]] && set -x
                
                set -o pipefail
                prev_total=0
                declare -a prev_time
                mem_tot=$(< /proc/meminfo grep MemTotal | awk '{print $2}')
                num_cpus=$(< /proc/cpuinfo grep '^processor' -c)
                
                nxf_pcpu() {
                    local pid=$1
                    local proc_time=$(2> /dev/null < /proc/$pid/stat awk '{sum=$14+$15; printf "%.0f",sum}' || echo -n 0)
                    local cpu_usage=$(echo -n $proc_time ${prev_time[pid]:-0} $total_time $prev_total $num_cpus | awk '{ pct=($1-$2)/($3-$4)*$5 *100; printf "%.1f", pct }' )
                    prev_time[pid]=$proc_time
                    nxf_pcpu_ret=$cpu_usage
                }
                
                nxf_tree() {
                    declare -a ALL_CHILD
                    while read P PP;do
                        ALL_CHILD[$PP]+=" $P"
                    done < <(ps -e -o pid= -o ppid=)
                
                    stat() {
                        nxf_pcpu $1 
                        local x_pid=$1
                        local x_stat=$(2> /dev/null < /proc/$1/stat awk '{print $3}' || echo -n X)
                        local x_pcpu=$nxf_pcpu_ret

                        local x_vsz=$(2> /dev/null < /proc/$1/stat awk '{printf "%.0f", $23/1024}' || echo -n 0)
                        local x_rss=$(2> /dev/null < /proc/$1/status grep VmRSS | awk '{print $2}' || echo -n 0)
                        local x_pmem=$(echo $x_rss | awk -v mem_tot=$mem_tot '{printf "%.1f", $1/mem_tot*100}')
                
                        local x_io=$(2> /dev/null < /proc/$1/io sed 's/^.*:\\s*//' | tr '\\n' ' ' || echo -n 0)
                        local x_vm=$(2> /dev/null < /proc/$1/status egrep 'VmPeak|VmHWM' | sed 's/^.*:\\s*//' | sed 's/[\\sa-zA-Z]*$//' | tr '\\n' ' ' || echo -n 0)
                
                        stat_ret+="$x_pid $x_stat $x_pcpu $x_pmem $x_vsz $x_rss"
                        if [[ $x_vm ]]; then stat_ret+=" $x_vm"; else stat_ret+=" 0 0"; fi
                        if [[ $x_io ]]; then stat_ret+=" $x_io"; fi
                        stat_ret+='\\n\'
                    }
                
                    walk() {
                        stat $1 
                        for i in ${ALL_CHILD[$1]:=}; do walk $i; done
                    }
                
                    stat_ret=\'\'
                    total_time=$(grep '^cpu ' /proc/stat |awk '{sum=$2+$3+$4+$5+$6+$7+$8+$9+$10; printf "%.0f",sum}')
                    walk $1
                    prev_total=$total_time
                    nxf_tree_ret=$stat_ret  
                }
                
                nxf_pstat() {
                    nxf_tree $1
                    if [[ "$nxf_tree_ret" ]]; then
                      nxf_pstat_ret=$(printf "$nxf_tree_ret" | awk '{ t3+=($3*10); t4+=($4*10); t5+=$5; t6+=$6; t7+=$7; t8+=$8; t9+=$9; t10+=$10; t11+=$11; t12+=$12; t13+=$13; t14+=$14 } END { printf "%d 0 %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f\\n", NR,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14 }')
                    else
                      nxf_pstat_ret=''  
                    fi
                }
                
                nxf_sleep() {
                  if [[ $1 < 0 ]]; then sleep 5;
                  elif [[ $1 < 10 ]]; then sleep 0.1 2>/dev/null || sleep 1;
                  elif [[ $1 < 130 ]]; then sleep 1;
                  else sleep 5; fi
                }
                
                nxf_date() {
                    local ts=$(date +%s%3N); [[ $ts == *3N ]] && date +%s000 || echo $ts
                }
                
                nxf_trace() {
                  local pid=$1; local trg=$2;
                  local count=0;
                  declare -a max=(); for i in {0..13}; do max[i]=0; done
                  while [[ -d /proc/$pid ]]; do
                    nxf_pstat $pid
                    if [[ "$nxf_pstat_ret" ]]; then
                    IFS=' ' read -a val <<< "$nxf_pstat_ret"; unset IFS
                    for i in {0..13}; do
                      [ ${val[i]} -gt ${max[i]} ] && max[i]=${val[i]}
                    done
                    echo "pid state %cpu %mem vmem rss peak_vmem peak_rss rchar wchar syscr syscw read_bytes write_bytes" > $trg
                    echo "${max[@]}" >> $trg
                    fi
                    nxf_sleep $count
                    count=$((count+1))
                  done
                } 
                
                
                trap 'exit ${ret:=$?}' EXIT
                touch .command.trace
                start_millis=$(nxf_date)
                (
                /bin/bash my-command.sh
                ) &
                pid=$!
                nxf_trace "$pid" .command.trace &
                mon=$!
                wait $pid || ret=$?
                end_millis=$(nxf_date)
                wait $mon
                echo $((end_millis-start_millis)) >> .command.trace
                
                # patch root ownership problem of files created with docker
                [ ${NXF_OWNER:=''} ] && chown -fR --from root $NXF_OWNER null/{*,.*} || true
                '''
                .stripIndent().leftTrim()

        cleanup:
        stub?.delete()
    }

    def 'should create container env' () {
        given:
        def bash = Spy(BashWrapperBuilder)

        when:
        def builder = bash.createContainerBuilder(null)
        then:
        bash.getEnvironment() >> [:]
        bash.getBinDir() >> Paths.get('/my/bin')
        bash.getWorkDir() >> Paths.get('/my/work/dir')
        bash.getStatsEnabled() >> false

        bash.getResolvedInputs() >> [:]
        bash.getContainerConfig() >> [engine: 'singularity', envWhitelist: 'FOO,BAR']
        bash.getContainerImage() >> 'foo/bar'
        bash.getContainerExecutable() >> false
        bash.getContainerMount() >> null
        bash.getContainerMemory() >> null
        bash.getContainerCpuset() >> null
        bash.getContainerOptions() >> null

        builder instanceof SingularityBuilder
        builder.env == ['FOO','BAR']
        builder.workDir == Paths.get('/my/work/dir')
    }

}
