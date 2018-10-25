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
        copy.touchFile(RUN) == "echo start | aws s3 cp --only-show-errors - s3://some/data/.command.run"
        copy.copyFile("nobel_prize_results.gz",Paths.get("/some/data/nobel_prize_results.gz")) == "nxf_s3_upload nobel_prize_results.gz s3://some/data"
        copy.exitFile(EXIT) == "| aws s3 cp --only-show-errors - s3://some/path/.exitcode || true"
        copy.stageInputFile(FILE, 'foo.txt') == """
                                    downloads+=("nxf_s3_download s3://some/data/nobel_prize_results.gz foo.txt")
                                    """
                                    .stripIndent().trim()
        copy.getUnstageOutputFilesScript(OUTPUTS,TARGET) == """
                                        uploads=()
                                        uploads+=("nxf_s3_upload 'outputs_*' s3://data/results")
                                        uploads+=("nxf_s3_upload 'final_folder' s3://data/results")
                                        nxf_parallel "\${uploads[@]}"
                                        """
                                        .stripIndent().trim()
    }

    def 'should return unstage script' () {
        given:
        def copy = new AwsBatchFileCopyStrategy(Mock(TaskBean), new AwsOptions())
        def target = Paths.get('/foo/bar')

        when:
        def script = copy.getUnstageOutputFilesScript(['file.txt'],target)
        then:
        script.trim() == """
                    uploads=()
                    uploads+=("nxf_s3_upload 'file.txt' s3://foo/bar")
                    nxf_parallel "\${uploads[@]}"
                    """
                    .stripIndent().trim()

        when:
        script = copy.getUnstageOutputFilesScript(['file-*.txt'],target)
        then:
        script.trim() == """
                        uploads=()
                        uploads+=("nxf_s3_upload 'file-*.txt' s3://foo/bar")
                        nxf_parallel "\${uploads[@]}"
                        """
                        .stripIndent().trim()

        when:
        script = copy.getUnstageOutputFilesScript(['file-[a,b].txt'],target)
        then:
        script.trim() == """
                    uploads=()
                    uploads+=("nxf_s3_upload 'file-[a,b].txt' s3://foo/bar")
                    nxf_parallel "\${uploads[@]}"
                    """
                    .stripIndent().trim()

        when:
        script = copy.getUnstageOutputFilesScript(['file-01(A).txt'],target)
        then:
        script.trim() == """
                    uploads=()
                    uploads+=("nxf_s3_upload 'file-01\\(A\\).txt' s3://foo/bar")
                    nxf_parallel "\${uploads[@]}"
                    """
                    .stripIndent().trim()
    }

    def 'should check the beforeScript' () {

        given:
        def bean = Mock(TaskBean)
        def opts = Mock(AwsOptions)
        def copy = Spy(AwsBatchFileCopyStrategy, constructorArgs: [bean, opts])

        when:
        def script = copy.getBeforeStartScript()
        then:
        1 * opts.getAwsCli() >> 'aws'
        1 * opts.getStorageClass() >> null
        1 * opts.getStorageEncryption() >> null

        script ==   '''
                    # aws helper
                    nxf_s3_upload() {
                        local pattern=$1
                        local s3path=$2
                        IFS=$'\\n'
                        for name in $(eval "ls -1d $pattern");do
                          if [[ -d "$name" ]]; then
                            aws s3 cp --only-show-errors --recursive --storage-class STANDARD "$name" "$s3path/$name"
                          else
                            aws s3 cp --only-show-errors --storage-class STANDARD "$name" "$s3path/$name"
                          fi
                        done
                        unset IFS
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
                    }     
                    '''.stripIndent()

        when:
        script = copy.getBeforeStartScript()
        then:
        1 * opts.getAwsCli() >> '/foo/aws'
        1 * opts.getStorageClass() >> 'REDUCED_REDUNDANCY'
        2 * opts.getStorageEncryption() >> 'AES256'

        script == '''
                # aws helper
                nxf_s3_upload() {
                    local pattern=$1
                    local s3path=$2
                    IFS=$'\\n'
                    for name in $(eval "ls -1d $pattern");do
                      if [[ -d "$name" ]]; then
                        /foo/aws s3 cp --only-show-errors --recursive --sse AES256 --storage-class REDUCED_REDUNDANCY "$name" "$s3path/$name"
                      else
                        /foo/aws s3 cp --only-show-errors --sse AES256 --storage-class REDUCED_REDUNDANCY "$name" "$s3path/$name"
                      fi
                    done
                    unset IFS
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
                
                nxf_parallel() {
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
                }     
            '''.stripIndent()
    }

    def 'should return env variables' () {

        given:
        def ENV = [FOO: 'hola', BAR:'world', PATH:'xxx']
        def bean = Mock(TaskBean)
        def opts = Mock(AwsOptions)
        def copy = Spy(AwsBatchFileCopyStrategy, constructorArgs: [bean, opts])

        when:
        def script = copy.getEnvScript(ENV)
        then:
        // note: PATH is always removed
        opts.getRemoteBinDir() >> null
        opts.getCliPath() >> null
        script == '''
            export BAR="world"
            export FOO="hola"
            '''.stripIndent().leftTrim()

        when:
        script = copy.getEnvScript(ENV)
        then:
        opts.getRemoteBinDir() >> '/foo/bar'
        opts.getAwsCli() >> 'aws'
        script == '''
            aws s3 cp --recursive --only-show-errors s3://foo/bar $PWD/nextflow-bin
            chmod +x $PWD/nextflow-bin/*
            export PATH=$PWD/nextflow-bin:$PATH
            export BAR="world"
            export FOO="hola"
            '''.stripIndent().leftTrim()

        when:
        script = copy.getEnvScript(ENV)
        then:
        opts.getAwsCli() >> '/conda/bin/aws'
        opts.getRemoteBinDir() >> '/foo/bar'
        script == '''
            /conda/bin/aws s3 cp --recursive --only-show-errors s3://foo/bar $PWD/nextflow-bin
            chmod +x $PWD/nextflow-bin/*
            export PATH=$PWD/nextflow-bin:$PATH
            export BAR="world"
            export FOO="hola"
            '''.stripIndent().leftTrim()

        when:
        script = copy.getEnvScript(ENV)
        then:
        opts.getAwsCli() >> '/conda/bin/aws'
        opts.getRemoteBinDir() >> '/foo/bar'
        opts.getRegion() >> 'eu-west-1'
        script == '''
            /conda/bin/aws s3 cp --recursive --only-show-errors s3://foo/bar $PWD/nextflow-bin
            chmod +x $PWD/nextflow-bin/*
            export PATH=$PWD/nextflow-bin:$PATH
            export BAR="world"
            export FOO="hola"
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
