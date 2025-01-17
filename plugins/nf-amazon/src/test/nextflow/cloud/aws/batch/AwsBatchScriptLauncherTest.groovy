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

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Paths

import nextflow.Session
import nextflow.SysEnv
import nextflow.cloud.aws.config.AwsConfig
import nextflow.cloud.aws.util.S3PathFactory
import nextflow.processor.TaskBean
import nextflow.util.Duration
import spock.lang.Specification
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
        def cfg = new AwsConfig(region: 'eu-west-1', batch: [cliPath:'/conda/bin/aws', retryMode: 'built-in'])
        def opts = new AwsOptions(awsConfig: cfg)
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
        binding.unstage_outputs == null

        binding.helpers_script == '''\
                # bash helper functions
                nxf_cp_retry() {
                    local max_attempts=5
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
                      /conda/bin/aws --region eu-west-1 s3 cp --only-show-errors --storage-class STANDARD - "$s3path"
                    elif [[ -d "$name" ]]; then
                      /conda/bin/aws --region eu-west-1 s3 cp --only-show-errors --recursive --storage-class STANDARD "$name" "$s3path/$name"
                    else
                      /conda/bin/aws --region eu-west-1 s3 cp --only-show-errors --storage-class STANDARD "$name" "$s3path/$name"
                    fi
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
                
                '''.stripIndent(true)
    }

    def 'should create task environment' () {
        /*
         * simple bash run
         */
        when:
        def bucket = Paths.get('/bucket/work')
        def opts = new AwsOptions(remoteBinDir: '/bucket/bin', awsConfig: new AwsConfig([:]))

        def binding = new AwsBatchScriptLauncher([
                name: 'Hello 1',
                workDir: bucket,
                targetDir: bucket,
                environment: [PATH:'/this:/that', FOO: 'xxx'],
                script: 'echo Hello world!' ] as TaskBean, opts) .makeBinding()

        then:
        binding.task_env == '''\
                    aws s3 cp --recursive --only-show-errors s3://bucket/bin $PWD/nextflow-bin
                    chmod +x $PWD/nextflow-bin/* || true
                    export PATH=$PWD/nextflow-bin:$PATH
                    export FOO="xxx"
                    '''.stripIndent()
    }

    def 'should cleanup temp files' () {

        when:
        def bucket = Paths.get('/bucket/work')
        def opts = new AwsOptions(remoteBinDir: '/bucket/bin', awsConfig: new AwsConfig([:]))

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
        def opts = new AwsOptions(awsConfig: new AwsConfig(batch: [retryMode: 'built-in']))

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
                downloads=(true)
                rm -f .command.sh
                rm -f .command.run
                rm -f .command.in
                downloads+=("nxf_cp_retry nxf_s3_download s3://bucket/work/.command.sh .command.sh")
                downloads+=("nxf_cp_retry nxf_s3_download s3://bucket/work/.command.run .command.run")
                downloads+=("nxf_cp_retry nxf_s3_download s3://bucket/work/.command.in .command.in")
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()

        binding.unstage_outputs == '''
                    uploads=()
                    IFS=$'\\n'
                    for name in $(eval "ls -1d foo.txt bar.fastq" | sort | uniq); do
                        uploads+=("nxf_s3_upload '$name' s3://bucket/work")
                    done
                    unset IFS
                    nxf_parallel "${uploads[@]}"
                    '''.stripIndent().leftTrim()

        binding.launch_cmd == '/bin/bash .command.run nxf_trace'
        
        binding.task_env == ''

        binding.helpers_script == '''\
                    # bash helper functions
                    nxf_cp_retry() {
                        local max_attempts=5
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


    def 'test bash wrapper with custom scratch'() {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def opts = new AwsOptions(awsConfig: new AwsConfig(aws:[batch:[cliPath:'/conda/bin/aws', region: 'eu-west-1']]))
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

    def 'test should disable scratch'() {

        given:
        def folder = Files.createTempDirectory('test')

        /*
         * simple bash run
         */
        when:
        def cfg = new AwsConfig(batch: [cliPath:'/conda/bin/aws'], region: 'eu-west-1')
        def opts = new AwsOptions(awsConfig: cfg)
        def bash = new AwsBatchScriptLauncher([
                name: 'Hello 1',
                workDir: folder,
                script: 'echo Hello world!',
                scratch: false
        ] as TaskBean, opts)
        bash.build()

        then:
        Files.exists(folder.resolve('.command.sh'))
        Files.exists(folder.resolve('.command.run'))

        folder.resolve('.command.run').text.contains("NXF_SCRATCH=''")

        cleanup:
        folder?.deleteDir()
    }

    def 'test download retry enabled'() {

        /*
         * simple bash run
         */
        when:
        def bucket = Paths.get('/bucket/work')
        def cfg = new AwsConfig(batch: [maxTransferAttempts:3, delayBetweenAttempts: '9 sec' as Duration, retryMode: 'built-in'])
        def opts = new AwsOptions(awsConfig: cfg)

        def binding = new AwsBatchScriptLauncher([
                name: 'Hello 1',
                workDir: bucket,
                // targetDir: bucket,
                script: 'echo Hello world!',
        ] as TaskBean, opts) .makeBinding()

        then:

        binding.stage_inputs == '''\
                # stage input files
                downloads=(true)
                rm -f .command.sh
                downloads+=("nxf_cp_retry nxf_s3_download s3://bucket/work/.command.sh .command.sh")
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()

        binding.helpers_script == '''\
                    # bash helper functions
                    nxf_cp_retry() {
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

    def 'should aws cli native retry'() {
        /*
         * simple bash run
         */
        when:
        def bucket = Paths.get('/bucket/work')
        def cfg = new AwsConfig(batch: [maxTransferAttempts: 3, retryMode: 'adaptive', delayBetweenAttempts: '9 sec' as Duration])
        def opts = new AwsOptions(awsConfig: cfg)

        def binding = new AwsBatchScriptLauncher([
                name: 'Hello 1',
                workDir: bucket,
                // targetDir: bucket,
                script: 'echo Hello world!',
        ] as TaskBean, opts) .makeBinding()

        then:

        binding.stage_inputs == '''\
                # stage input files
                downloads=(true)
                rm -f .command.sh
                downloads+=("nxf_s3_download s3://bucket/work/.command.sh .command.sh")
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()

        binding.helpers_script == '''\
                    # bash helper functions
                    nxf_cp_retry() {
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
                    
                    # aws cli retry config
                    export AWS_RETRY_MODE=adaptive 
                    export AWS_MAX_ATTEMPTS=3
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


    def 'should include fix ownership command' () {
        given:
        def cfg = new AwsConfig(batch: [cliPath:'/conda/bin/aws'], region: 'eu-west-1')
        def opts = new AwsOptions(awsConfig: cfg)
        def builder = new AwsBatchScriptLauncher([
                name: 'Hello 1',
                workDir: Paths.get('/work/dir'),
                script: 'echo Hello world!',
                containerConfig: [fixOwnership: true],
                input: 'Ciao ciao' ] as TaskBean, opts)

        when:
        def binding = builder.makeBinding()
        then:
        builder.fixOwnership() >> true
        binding.fix_ownership == '[ ${NXF_OWNER:=\'\'} ] && (shopt -s extglob; GLOBIGNORE=\'..\'; chown -fR --from root $NXF_OWNER /work/dir/{*,.*}) || true'

    }

    def 'should not create separate stage script' () {
        given:
        SysEnv.push([NXF_WRAPPER_STAGE_FILE_THRESHOLD: '100'])
        and:
        def workDir = S3PathFactory.parse('s3://my-bucket/work')
        and:
        def inputFiles = [
                'sample_1.fq': Paths.get('/my-bucket/data/sample_1.fq'),
                'sample_2.fq': Paths.get('/my-bucket/data/sample_2.fq'),
        ]
        def stageScript = '''\
                # stage input files
                downloads=(true)
                rm -f sample_1.fq
                rm -f sample_2.fq
                rm -f .command.sh
                downloads+=("nxf_s3_download s3://my-bucket/data/sample_1.fq sample_1.fq")
                downloads+=("nxf_s3_download s3://my-bucket/data/sample_2.fq sample_2.fq")
                downloads+=("nxf_s3_download s3://my-bucket/work/.command.sh .command.sh")
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()
        and:
        def bean = [
                workDir: workDir,
                targetDir: workDir,
                inputFiles: inputFiles,
                script: 'echo Hello world!'
        ] as TaskBean
        def opts = new AwsOptions()
        def builder = new AwsBatchScriptLauncher(bean, opts)

        when:
        def binding = builder.makeBinding()
        then:
        binding.stage_inputs == stageScript

        cleanup:
        SysEnv.pop()
    }

}
