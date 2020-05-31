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

package nextflow.cloud.aws.batch


import nextflow.util.Duration
import spock.lang.Specification

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.Session
import nextflow.processor.TaskBean
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsBatchScriptLauncherTest extends Specification {

    def setup() {
        new Session()
    }

    def 'test bash wrapper with input'() {

        /*
         * simple bash run
         */
        when:
        def opts = new AwsOptions(cliPath:'/conda/bin/aws', region: 'eu-west-1')
        def binding = new AwsBatchScriptLauncher([
                name: 'Hello 1',
                workDir: Paths.get('/work/dir'),
                script: 'echo Hello world!',
                environment: [FOO: 1, BAR:'any'],
                input: 'Ciao ciao' ] as TaskBean, opts) .makeBinding()

        then:
        binding.unstage_controls == '''\
                nxf_s3_upload .command.out s3://work/dir || true
                nxf_s3_upload .command.err s3://work/dir || true
                '''.stripIndent()

        binding.launch_cmd == '/bin/bash -ue .command.sh < .command.in'
        binding.unstage_outputs == ''

        binding.helpers_script == '''
                # aws helper
                nxf_s3_upload() {
                    local pattern=$1
                    local s3path=$2
                    IFS=$'\\n\'
                    for name in $(eval "ls -1d $pattern");do
                      if [[ -d "$name" ]]; then
                        /conda/bin/aws --region eu-west-1 s3 cp --only-show-errors --recursive --storage-class STANDARD "$name" "$s3path/$name"
                      else
                        /conda/bin/aws --region eu-west-1 s3 cp --only-show-errors --storage-class STANDARD "$name" "$s3path/$name"
                      fi
                    done
                    unset IFS
                }
                
                nxf_s3_retry() {
                    local max_attempts=1
                    local timeout=10
                    local attempt=0
                    local exitCode=0
                    while (( \$attempt < \$max_attempts ))
                    do
                      if "\$@"
                        then
                          return 0
                      else
                        exitCode=\$?
                      fi
                      if [[ \$exitCode == 0 ]]
                      then
                        break
                      fi
                      sleep \$timeout
                      attempt=\$(( attempt + 1 ))
                      timeout=\$(( timeout * 2 ))
                    done
                }
                
                nxf_s3_download() {
                    local source=$1
                    local target=$2
                    local file_name=$(basename $1)
                    local is_dir=$(/conda/bin/aws --region eu-west-1 s3 ls $source | grep -F "PRE ${file_name}/" -c)
                    if [[ $is_dir == 1 ]]; then
                        /conda/bin/aws --region eu-west-1 s3 cp --only-show-errors --recursive "$source" "$target"
                    else 
                        /conda/bin/aws --region eu-west-1 s3 cp --only-show-errors "$source" "$target"
                    fi
                }
                
                nxf_parallel() {
                    IFS=$'\\n\'
                    local cmd=("$@")
                    local cpus=$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
                    local max=$(if (( cpus>16 )); then echo 16; else echo $cpus; fi)
                    local i=0
                    local pid=()
                    (
                    set +u
                    while ((i<${#cmd[@]})); do
                        local copy=()
                        for x in "${pid[@]}"; do
                          [[ -e /proc/$x ]] && copy+=($x) 
                        done
                        pid=("${copy[@]}")
                
                        if ((${#pid[@]}>=$max)); then 
                          sleep 1 
                        else 
                          eval "${cmd[$i]}" &
                          pid+=($!)
                          ((i+=1))
                        fi 
                    done
                    ((${#pid[@]}>0)) && wait ${pid[@]}
                    )
                    unset IFS
                }
                
                '''.stripIndent()
    }

    def 'should create task environment' () {
        /*
         * simple bash run
         */
        when:
        def bucket = Paths.get('/bucket/work')
        def opts = new AwsOptions(remoteBinDir: '/bucket/bin')

        def binding = new AwsBatchScriptLauncher([
                name: 'Hello 1',
                workDir: bucket,
                targetDir: bucket,
                environment: [PATH:'/this:/that', FOO: 'xxx'],
                script: 'echo Hello world!' ] as TaskBean, opts) .makeBinding()

        then:
        binding.task_env == '''\
                    aws s3 cp --recursive --only-show-errors s3://bucket/bin $PWD/nextflow-bin
                    chmod +x $PWD/nextflow-bin/*
                    export PATH=$PWD/nextflow-bin:$PATH
                    export FOO="xxx"
                    '''.stripIndent()
    }

    def 'should cleanup temp files' () {

        when:
        def bucket = Paths.get('/bucket/work')
        def opts = new AwsOptions(remoteBinDir: '/bucket/bin')

        def binding = new AwsBatchScriptLauncher([
                name: 'Hello 1',
                workDir: bucket,
                targetDir: bucket,
                environment: [PATH:'/this:/that', FOO: 'xxx'],
                script: 'echo Hello world!' ] as TaskBean, opts) .makeBinding()

        then:
        binding.cleanup_cmd == 'rm -rf $NXF_SCRATCH || true\n'
    }

    def 'test bash wrapper with outputs and stats'() {

        /*
         * simple bash run
         */
        when:
        def bucket = Paths.get('/bucket/work')
        def opts = new AwsOptions()

        def binding = new AwsBatchScriptLauncher([
                name: 'Hello 1',
                workDir: bucket,
                targetDir: bucket,
                statsEnabled: true,
                outputFiles: ['foo.txt', 'bar.fastq'],
                script: 'echo Hello world!',
                input: 'Ciao ciao' ] as TaskBean, opts) .makeBinding()

        then:

        binding.unstage_controls == '''\
                nxf_s3_upload .command.out s3://bucket/work || true
                nxf_s3_upload .command.err s3://bucket/work || true
                nxf_s3_upload .command.trace s3://bucket/work || true
                '''.stripIndent()

        binding.stage_inputs == '''\
                # stage input files
                downloads=()
                rm -f .command.sh
                rm -f .command.run
                rm -f .command.in
                downloads+=("nxf_s3_download s3://bucket/work/.command.sh .command.sh")
                downloads+=("nxf_s3_download s3://bucket/work/.command.run .command.run")
                downloads+=("nxf_s3_download s3://bucket/work/.command.in .command.in")
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()

        binding.unstage_outputs == '''\
                  uploads=()
                  uploads+=("nxf_s3_upload 'foo.txt' s3://bucket/work")
                  uploads+=("nxf_s3_upload 'bar.fastq' s3://bucket/work")
                  nxf_parallel "${uploads[@]}"
                  '''.stripIndent().rightTrim()

        binding.launch_cmd == '/bin/bash .command.run nxf_trace'
        
        binding.task_env == ''

        binding.helpers_script == '''
                    # aws helper
                    nxf_s3_upload() {
                        local pattern=$1
                        local s3path=$2
                        IFS=$'\\n\'
                        for name in $(eval "ls -1d $pattern");do
                          if [[ -d "$name" ]]; then
                            aws s3 cp --only-show-errors --recursive --storage-class STANDARD "$name" "$s3path/$name"
                          else
                            aws s3 cp --only-show-errors --storage-class STANDARD "$name" "$s3path/$name"
                          fi
                        done
                        unset IFS
                    }
                    
                    nxf_s3_retry() {
                        local max_attempts=1
                        local timeout=10
                        local attempt=0
                        local exitCode=0
                        while (( \$attempt < \$max_attempts ))
                        do
                          if "\$@"
                            then
                              return 0
                          else
                            exitCode=\$?
                          fi
                          if [[ \$exitCode == 0 ]]
                          then
                            break
                          fi
                          sleep \$timeout
                          attempt=\$(( attempt + 1 ))
                          timeout=\$(( timeout * 2 ))
                        done
                    }
                    
                    nxf_s3_download() {
                        local source=$1
                        local target=$2
                        local file_name=$(basename $1)
                        local is_dir=$(aws s3 ls $source | grep -F "PRE ${file_name}/" -c)
                        if [[ $is_dir == 1 ]]; then
                            aws s3 cp --only-show-errors --recursive "$source" "$target"
                        else 
                            aws s3 cp --only-show-errors "$source" "$target"
                        fi
                    }
                    
                    nxf_parallel() {
                        IFS=$'\\n\'
                        local cmd=("$@")
                        local cpus=$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
                        local max=$(if (( cpus>16 )); then echo 16; else echo $cpus; fi)
                        local i=0
                        local pid=()
                        (
                        set +u
                        while ((i<${#cmd[@]})); do
                            local copy=()
                            for x in "${pid[@]}"; do
                              [[ -e /proc/$x ]] && copy+=($x) 
                            done
                            pid=("${copy[@]}")
                    
                            if ((${#pid[@]}>=$max)); then 
                              sleep 1 
                            else 
                              eval "${cmd[$i]}" &
                              pid+=($!)
                              ((i+=1))
                            fi 
                        done
                        ((${#pid[@]}>0)) && wait ${pid[@]}
                        )
                        unset IFS
                    }
                    
                    '''.stripIndent()

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

    def 'test download retry enabled'() {

        /*
         * simple bash run
         */
        when:
        def bucket = Paths.get('/bucket/work')
        def opts = new AwsOptions()
        opts.maxTransferAttempts = 3
        opts.delayBetweenAttempts = '9 sec' as Duration

        def binding = new AwsBatchScriptLauncher([
                name: 'Hello 1',
                workDir: bucket,
                // targetDir: bucket,
                script: 'echo Hello world!',
        ] as TaskBean, opts) .makeBinding()

        then:

        binding.stage_inputs == '''\
                # stage input files
                downloads=()
                rm -f .command.sh
                downloads+=("nxf_s3_retry nxf_s3_download s3://bucket/work/.command.sh .command.sh")
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()

        binding.helpers_script == '''
                    # aws helper
                    nxf_s3_upload() {
                        local pattern=$1
                        local s3path=$2
                        IFS=$'\\n\'
                        for name in $(eval "ls -1d $pattern");do
                          if [[ -d "$name" ]]; then
                            aws s3 cp --only-show-errors --recursive --storage-class STANDARD "$name" "$s3path/$name"
                          else
                            aws s3 cp --only-show-errors --storage-class STANDARD "$name" "$s3path/$name"
                          fi
                        done
                        unset IFS
                    }
                    
                    nxf_s3_retry() {
                        local max_attempts=3
                        local timeout=9
                        local attempt=0
                        local exitCode=0
                        while (( \$attempt < \$max_attempts ))
                        do
                          if "\$@"
                            then
                              return 0
                          else
                            exitCode=\$?
                          fi
                          if [[ \$exitCode == 0 ]]
                          then
                            break
                          fi
                          sleep \$timeout
                          attempt=\$(( attempt + 1 ))
                          timeout=\$(( timeout * 2 ))
                        done
                    }
                    
                    nxf_s3_download() {
                        local source=$1
                        local target=$2
                        local file_name=$(basename $1)
                        local is_dir=$(aws s3 ls $source | grep -F "PRE ${file_name}/" -c)
                        if [[ $is_dir == 1 ]]; then
                            aws s3 cp --only-show-errors --recursive "$source" "$target"
                        else 
                            aws s3 cp --only-show-errors "$source" "$target"
                        fi
                    }
                    
                    nxf_parallel() {
                        IFS=$'\\n\'
                        local cmd=("$@")
                        local cpus=$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
                        local max=$(if (( cpus>16 )); then echo 16; else echo $cpus; fi)
                        local i=0
                        local pid=()
                        (
                        set +u
                        while ((i<${#cmd[@]})); do
                            local copy=()
                            for x in "${pid[@]}"; do
                              [[ -e /proc/$x ]] && copy+=($x) 
                            done
                            pid=("${copy[@]}")
                    
                            if ((${#pid[@]}>=$max)); then 
                              sleep 1 
                            else 
                              eval "${cmd[$i]}" &
                              pid+=($!)
                              ((i+=1))
                            fi 
                        done
                        ((${#pid[@]}>0)) && wait ${pid[@]}
                        )
                        unset IFS
                    }
                    
                    '''.stripIndent()

    }

}
