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

package nextflow.cloud.aws.batch

import java.nio.file.Files

import nextflow.processor.TaskBean
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsBatchScriptLauncherTest extends Specification {

    def 'test bash wrapper with input'() {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def opts = new AwsOptions(cliPath:'/conda/bin/aws', region: 'eu-west-1')
        def bash = new AwsBatchScriptLauncher([
                name: 'Hello 1',
                workDir: folder,
                script: 'echo Hello world!',
                environment: [FOO: 1, BAR:'any'],
                input: 'Ciao ciao'
        ] as TaskBean, opts)
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
                  printf \$exit_status | /conda/bin/aws --region eu-west-1 s3 cp --only-show-errors - s3:/${folder}/.exitcode || true
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
                        /conda/bin/aws --region eu-west-1 s3 cp --only-show-errors --recursive --storage-class STANDARD "\$name" "\$s3path/\$name"
                      else
                        /conda/bin/aws --region eu-west-1 s3 cp --only-show-errors --storage-class STANDARD "\$name" "\$s3path/\$name"
                      fi
                    done
                    unset IFS
                }

                nxf_s3_download() {
                    local source=\$1
                    local target=\$2
                    local file_name=\$(basename \$1)
                    local is_dir=\$(/conda/bin/aws --region eu-west-1 s3 ls \$source | grep -F "PRE \${file_name}/" -c)
                    if [[ \$is_dir == 1 ]]; then
                        /conda/bin/aws --region eu-west-1 s3 cp --only-show-errors --recursive "\$source" "\$target"
                    else 
                        /conda/bin/aws --region eu-west-1 s3 cp --only-show-errors "\$source" "\$target"
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

                echo start | /conda/bin/aws --region eu-west-1 s3 cp --only-show-errors - s3:/${folder}/.command.begin
                # task environment
                export FOO="1"
                export BAR="any"

                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH
                # stage input files
                downloads=()
                rm -f .command.sh
                rm -f .command.in
                downloads+=("nxf_s3_download s3:/$folder/.command.sh .command.sh")
                downloads+=("nxf_s3_download s3:/$folder/.command.in .command.in")
                nxf_parallel "\${downloads[@]}"


                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                /bin/bash -ue .command.sh < .command.in
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                nxf_s3_upload .command.out s3:/${folder} || true
                nxf_s3_upload .command.err s3:/${folder} || true
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }


    def 'test bash wrapper with outputs and stats'() {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def opts = new AwsOptions()
        def bash = new AwsBatchScriptLauncher([
                name: 'Hello 1',
                workDir: folder,
                targetDir: folder,
                statsEnabled: true,
                outputFiles: ['foo.txt', 'bar.fastq'],
                script: 'echo Hello world!',
                input: 'Ciao ciao'
        ] as TaskBean, opts)
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
                  printf \$exit_status | aws s3 cp --only-show-errors - s3:/${folder}/.exitcode || true
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

                echo start | aws s3 cp --only-show-errors - s3:/${folder}/.command.begin
                [[ \$NXF_SCRATCH ]] && echo "nxf-scratch-dir \$HOSTNAME:\$NXF_SCRATCH" && cd \$NXF_SCRATCH
                # stage input files
                downloads=()
                rm -f .command.sh
                rm -f .command.stub
                rm -f .command.in
                downloads+=("nxf_s3_download s3:/$folder/.command.sh .command.sh")
                downloads+=("nxf_s3_download s3:/$folder/.command.stub .command.stub")
                downloads+=("nxf_s3_download s3:/$folder/.command.in .command.in")
                nxf_parallel "\${downloads[@]}"


                set +e
                ctmp=\$(set +u; nxf_mktemp /dev/shm 2>/dev/null || nxf_mktemp \$TMPDIR)
                cout=\$ctmp/.command.out; mkfifo \$cout
                cerr=\$ctmp/.command.err; mkfifo \$cerr
                tee .command.out < \$cout &
                tee1=\$!
                tee .command.err < \$cerr >&2 &
                tee2=\$!
                (
                /bin/bash .command.stub
                ) >\$cout 2>\$cerr &
                pid=\$!
                wait \$pid || ret=\$?
                wait \$tee1 \$tee2
                nxf_s3_upload .command.out s3:/${folder} || true
                nxf_s3_upload .command.err s3:/${folder} || true
                # copies output files to target
                if [[ \${ret:=0} == 0 ]]; then
                  uploads=()
                  uploads+=("nxf_s3_upload 'foo.txt' s3:/$folder")
                  uploads+=("nxf_s3_upload 'bar.fastq' s3:/$folder")
                  nxf_parallel "\${uploads[@]}"
                fi
                nxf_s3_upload .command.trace s3:/${folder} || true
                """
                        .stripIndent().leftTrim()


        cleanup:
        folder?.deleteDir()
    }


    def 'test bash wrapper with custom scratch'() {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def opts = new AwsOptions(cliPath:'/conda/bin/aws', region: 'eu-west-1')
        def bash = new AwsBatchScriptLauncher([
                name: 'Hello 1',
                workDir: folder,
                script: 'echo Hello world!',
                scratch: '/foo/bar/tmp'
        ] as TaskBean, opts)
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))

        folder.resolve('.command.run').text.contains('NXF_SCRATCH="$(set +u; nxf_mktemp /foo/bar/tmp)"')


        cleanup:
        folder?.deleteDir()
    }

}
