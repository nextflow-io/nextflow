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

package nextflow.executor

import spock.lang.Specification

import java.nio.file.Path
import java.nio.file.Paths

import nextflow.Global
import nextflow.Session
import nextflow.processor.TaskBean
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BashWrapperBuilderS3Test extends Specification {

    def 'should include s3 helpers' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [:] }
        and:
        def folder = Paths.get('/work/dir')
        def target = Mock(Path)
        target.toString() >> '/some/bucket'

        def bean = new TaskBean([
                name: 'Hello 1',
                workDir: folder,
                targetDir: target,
                scratch: true,
                outputFiles: ['test.bam','test.bai'],
                script: 'echo Hello world!',
        ])

        SimpleFileCopyStrategy copy = Spy(SimpleFileCopyStrategy, constructorArgs:[bean])
        copy.getPathScheme(target) >> 's3'

        /*
         * simple bash run
         */
        when:
        def binding = new BashWrapperBuilder(bean,copy).makeBinding()
        then:
        binding.unstage_outputs == '''\
                  nxf_s3_upload 'test.bam' s3://some/bucket || true
                  nxf_s3_upload 'test.bai' s3://some/bucket || true
                  '''.stripIndent().rightTrim()

        binding.helpers_script == '''\
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


}
