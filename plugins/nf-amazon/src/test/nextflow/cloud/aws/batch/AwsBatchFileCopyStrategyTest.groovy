/*
 * Copyright 2013-2024, Seqera Labs
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
import java.nio.file.Paths

import nextflow.processor.TaskBean
import spock.lang.Specification
import test.TestHelper

class AwsBatchFileCopyStrategyTest extends Specification {

    def 'should strip out file/folder name from target S3 path' () {
        given:
        def OUTPUTS = ["outputs_*","final_folder"]
        def TARGET =  Paths.get('/data/results')
        def FILE = Paths.get('/some/data/nobel_prize_results.gz')
        def EXIT = Paths.get('/some/path/.exitcode')
        def RUN = Paths.get('/some/data/.command.run')
        def copy = new AwsBatchFileCopyStrategy(Mock(TaskBean), new AwsOptions())
        expect:
        copy.touchFile(RUN) == "echo start | nxf_s3_upload - s3://some/data/.command.run"
        copy.copyFile("nobel_prize_results.gz",Paths.get("/some/data/nobel_prize_results.gz")) == "nxf_s3_upload nobel_prize_results.gz s3://some/data"
        copy.exitFile(EXIT) == "| nxf_s3_upload - s3://some/path/.exitcode || true"
        copy.stageInputFile(FILE, 'foo.txt') == """
                                    downloads+=("nxf_s3_download s3://some/data/nobel_prize_results.gz foo.txt")
                                    """
                                    .stripIndent().trim()
        copy.getUnstageOutputFilesScript(OUTPUTS,TARGET) == '''
                                        uploads=()
                                        IFS=$'\\n'
                                        for name in $(eval "ls -1d outputs_* final_folder" | sort | uniq); do
                                            uploads+=("nxf_s3_upload '$name' s3://data/results")
                                        done
                                        unset IFS
                                        nxf_parallel "${uploads[@]}"
                                        '''
                                        .stripIndent().leftTrim()
    }

    def 'should return unstage script' () {
        given:
        def copy = new AwsBatchFileCopyStrategy(Mock(TaskBean), new AwsOptions())
        def target = Paths.get('/foo/bar')

        when:
        def script = copy.getUnstageOutputFilesScript(['file.txt'],target)
        then:
        script.trim() == '''
                    uploads=()
                    IFS=$'\\n'
                    for name in $(eval "ls -1d file.txt" | sort | uniq); do
                        uploads+=("nxf_s3_upload '$name' s3://foo/bar")
                    done
                    unset IFS
                    nxf_parallel "${uploads[@]}"
                    '''
                    .stripIndent().trim()

        when:
        script = copy.getUnstageOutputFilesScript(['file-*.txt'],target)
        then:
        script.trim() == '''
                        uploads=()
                        IFS=$'\\n'
                        for name in $(eval "ls -1d file-*.txt" | sort | uniq); do
                            uploads+=("nxf_s3_upload '$name' s3://foo/bar")
                        done
                        unset IFS
                        nxf_parallel "${uploads[@]}"
                        '''
                        .stripIndent().trim()

        when:
        script = copy.getUnstageOutputFilesScript(['file-[a,b].txt'],target)
        then:
        script.trim() == '''
                    uploads=()
                    IFS=$'\\n'
                    for name in $(eval "ls -1d file-[a,b].txt" | sort | uniq); do
                        uploads+=("nxf_s3_upload '$name' s3://foo/bar")
                    done
                    unset IFS
                    nxf_parallel "${uploads[@]}"
                    '''
                    .stripIndent().trim()

        when:
        script = copy.getUnstageOutputFilesScript(['file-01(A).txt', 'f o o.txt'],target)
        then:
        script.trim() == '''
                    uploads=()
                    IFS=$'\\n'
                    for name in $(eval "ls -1d file-01\\(A\\).txt f\\ o\\ o.txt" | sort | uniq); do
                        uploads+=("nxf_s3_upload '$name' s3://foo/bar")
                    done
                    unset IFS
                    nxf_parallel "${uploads[@]}"
                    '''
                    .stripIndent().trim()
    }

    def 'should check the beforeScript' () {

        given:
        def bean = Mock(TaskBean)
        def opts = Mock(AwsOptions)
        AwsBatchFileCopyStrategy copy = Spy(AwsBatchFileCopyStrategy, constructorArgs: [bean, opts])

        when:
        def script = copy.getBeforeStartScript()
        then:
        1 * opts.getAwsCli() >> 'aws'
        1 * opts.getStorageClass() >> null
        1 * opts.getStorageEncryption() >> null

        script ==   '''\
                    # bash helper functions
                    nxf_cp_retry() {
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
                          nxf_sleep \$timeout
                          attempt=\$(( attempt + 1 ))
                          timeout=\$(( timeout * 2 ))
                        done
                    }
                    
                    nxf_parallel() {
                        IFS=$'\\n'
                        local cmd=("$@")
                        local cpus=$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
                        local max=$(if (( cpus>4 )); then echo 4; else echo $cpus; fi)
                        local i=0
                        local pid=()
                        (
                        set +u
                        while ((i<${#cmd[@]})); do
                            local copy=()
                            for x in "${pid[@]}"; do
                              # if the process exist, keep in the 'copy' array, otherwise wait on it to capture the exit code
                              # see https://github.com/nextflow-io/nextflow/pull/4050
                              [[ -e /proc/$x ]] && copy+=($x) || wait $x
                            done
                            pid=("${copy[@]}")
                    
                            if ((${#pid[@]}>=$max)); then
                              nxf_sleep 0.2
                            else
                              eval "${cmd[$i]}" &
                              pid+=($!)
                              ((i+=1))
                            fi
                        done
                        for p in "${pid[@]}"; do
                            wait $p
                        done
                        )
                        unset IFS
                    }
                    
                    # aws helper
                    nxf_s3_upload() {
                        local name=$1
                        local s3path=$2
                        if [[ "$name" == - ]]; then
                          aws s3 cp --only-show-errors --storage-class STANDARD - "$s3path"
                        elif [[ -d "$name" ]]; then
                          aws s3 cp --only-show-errors --recursive --storage-class STANDARD "$name" "$s3path/$name"
                        else
                          aws s3 cp --only-show-errors --storage-class STANDARD "$name" "$s3path/$name"
                        fi
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
                    '''.stripIndent(true)

        when:
        script = copy.getBeforeStartScript()
        then:
        1 * opts.getAwsCli() >> '/foo/aws'
        1 * opts.getStorageClass() >> 'STANDARD_IA'
        1 * opts.getStorageEncryption() >> 'AES256'

        script == '''\
                # bash helper functions
                nxf_cp_retry() {
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
                      nxf_sleep \$timeout
                      attempt=\$(( attempt + 1 ))
                      timeout=\$(( timeout * 2 ))
                    done
                }
                
                nxf_parallel() {
                    IFS=$'\\n'
                    local cmd=("$@")
                    local cpus=$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
                    local max=$(if (( cpus>4 )); then echo 4; else echo $cpus; fi)
                    local i=0
                    local pid=()
                    (
                    set +u
                    while ((i<${#cmd[@]})); do
                        local copy=()
                        for x in "${pid[@]}"; do
                          # if the process exist, keep in the 'copy' array, otherwise wait on it to capture the exit code
                          # see https://github.com/nextflow-io/nextflow/pull/4050
                          [[ -e /proc/$x ]] && copy+=($x) || wait $x
                        done
                        pid=("${copy[@]}")
                
                        if ((${#pid[@]}>=$max)); then
                          nxf_sleep 0.2
                        else
                          eval "${cmd[$i]}" &
                          pid+=($!)
                          ((i+=1))
                        fi
                    done
                    for p in "${pid[@]}"; do
                        wait $p
                    done
                    )
                    unset IFS
                }
                
                # aws helper
                nxf_s3_upload() {
                    local name=$1
                    local s3path=$2
                    if [[ "$name" == - ]]; then
                      /foo/aws s3 cp --only-show-errors --sse AES256 --storage-class STANDARD_IA - "$s3path"
                    elif [[ -d "$name" ]]; then
                      /foo/aws s3 cp --only-show-errors --recursive --sse AES256 --storage-class STANDARD_IA "$name" "$s3path/$name"
                    else
                      /foo/aws s3 cp --only-show-errors --sse AES256 --storage-class STANDARD_IA "$name" "$s3path/$name"
                    fi
                }
 
                nxf_s3_download() {
                    local source=$1
                    local target=$2
                    local file_name=$(basename $1)
                    local is_dir=$(/foo/aws s3 ls $source | grep -F "PRE ${file_name}/" -c)
                    if [[ $is_dir == 1 ]]; then
                        /foo/aws s3 cp --only-show-errors --recursive "$source" "$target"
                    else 
                        /foo/aws s3 cp --only-show-errors "$source" "$target"
                    fi
                }
            '''.stripIndent(true)
    }

    def 'should return env variables' () {

        given:
        def ENV = [FOO: 'hola', BAR:'world', PATH:'xxx']
        def bean = Mock(TaskBean)
        def opts = Mock(AwsOptions)
        AwsBatchFileCopyStrategy copy = Spy(AwsBatchFileCopyStrategy, constructorArgs: [bean, opts])

        when:
        def script = copy.getEnvScript(ENV,false)
        then:
        // note: PATH is always removed
        opts.getRemoteBinDir() >> null
        opts.getCliPath() >> null
        script == '''
            export FOO="hola"
            export BAR="world"
            '''.stripIndent().leftTrim()

        when:
        script = copy.getEnvScript(ENV,false)
        then:
        opts.getRemoteBinDir() >> '/foo/bar'
        opts.getAwsCli() >> 'aws'
        script == '''
            aws s3 cp --recursive --only-show-errors s3://foo/bar $PWD/nextflow-bin
            chmod +x $PWD/nextflow-bin/* || true
            export PATH=$PWD/nextflow-bin:$PATH
            export FOO="hola"
            export BAR="world"
            '''.stripIndent().leftTrim()

        when:
        script = copy.getEnvScript(ENV,false)
        then:
        opts.getAwsCli() >> '/conda/bin/aws'
        opts.getRemoteBinDir() >> '/foo/bar'
        script == '''
            /conda/bin/aws s3 cp --recursive --only-show-errors s3://foo/bar $PWD/nextflow-bin
            chmod +x $PWD/nextflow-bin/* || true
            export PATH=$PWD/nextflow-bin:$PATH
            export FOO="hola"
            export BAR="world"
            '''.stripIndent().leftTrim()

        when:
        script = copy.getEnvScript(ENV,false)
        then:
        opts.getAwsCli() >> '/conda/bin/aws'
        opts.getRemoteBinDir() >> '/foo/bar'
        opts.getRegion() >> 'eu-west-1'
        script == '''
            /conda/bin/aws s3 cp --recursive --only-show-errors s3://foo/bar $PWD/nextflow-bin
            chmod +x $PWD/nextflow-bin/* || true
            export PATH=$PWD/nextflow-bin:$PATH
            export FOO="hola"
            export BAR="world"
            '''.stripIndent().leftTrim()

    }


    def 'should return stage input input file'() {
        given:
        def file = TestHelper.createInMemTempFile('foo.txt')

        def bean = Mock(TaskBean)
        def opts = Mock(AwsOptions)
        def copy = new AwsBatchFileCopyStrategy(bean, opts)

        when:
        def script = copy.stageInputFile( file, 'bar.txt')
        then:
        script == "downloads+=(\"nxf_s3_download s3:/$file bar.txt\")" as String

    }
    
}
