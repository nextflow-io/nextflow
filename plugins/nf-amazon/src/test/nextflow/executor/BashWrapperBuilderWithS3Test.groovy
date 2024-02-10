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

package nextflow.executor

import java.nio.file.Paths

import nextflow.Global
import nextflow.Session
import nextflow.cloud.aws.batch.AwsBatchFileCopyStrategy
import nextflow.cloud.aws.batch.AwsOptions
import nextflow.cloud.aws.util.S3PathFactory
import nextflow.processor.TaskBean
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BashWrapperBuilderWithS3Test extends Specification {

    def 'should include s3 helpers' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [:] }
        and:
        def folder = Paths.get('/work/dir')
        def target = S3PathFactory.parse('s3://some/buck et')  // <-- path with blank

        def bean = new TaskBean([
                name: 'Hello 1',
                workDir: folder,
                targetDir: target,
                scratch: true,
                outputFiles: ['test.bam','test.bai', 'bla nk.txt'],  // <-- file name with blank
                script: 'echo Hello world!',
        ])

        def copy = new SimpleFileCopyStrategy(bean)

        /*
         * simple bash run
         */
        when:
        def binding = new BashWrapperBuilder(bean,copy).makeBinding()
        then:
        binding.unstage_outputs == '''\
                    IFS=$'\\n'
                    for name in $(eval "ls -1d test.bam test.bai bla\\ nk.txt" | sort | uniq); do
                        nxf_s3_upload $name s3://some/buck\\ et || true
                    done
                    unset IFS
                    '''.stripIndent().rightTrim()

        binding.helpers_script == '''\
            # aws cli retry config
            export AWS_RETRY_MODE=standard 
            export AWS_MAX_ATTEMPTS=5
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
    }

    def 'should include s3 helpers and bash lib' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [:] }
        and:
        def folder = Paths.get('/work/dir')
        def target = S3PathFactory.parse('s3://some/bucket')

        def bean = new TaskBean([
                name: 'Hello 1',
                workDir: folder,
                targetDir: target,
                scratch: true,
                outputFiles: ['test.bam','test.bai'],
                script: 'echo Hello world!',
        ])

        def copy = new AwsBatchFileCopyStrategy(bean, Mock(AwsOptions))

        /*
         * simple bash run
         */
        when:
        def binding = new BashWrapperBuilder(bean,copy).makeBinding()
        then:
        binding.unstage_outputs == '''\
                    uploads=()
                    IFS=$'\\n'
                    for name in $(eval "ls -1d test.bam test.bai" | sort | uniq); do
                        uploads+=("nxf_s3_upload '$name' s3://some/bucket")
                    done
                    unset IFS
                    nxf_parallel "${uploads[@]}"
                    '''.stripIndent()

        binding.helpers_script == '''\
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
                IFS=$'\\n\'
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
    }
}
