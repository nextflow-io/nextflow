package nextflow.cloud.aws.util

import nextflow.Global
import nextflow.Session
import nextflow.cloud.aws.batch.AwsOptions
import nextflow.util.Duration
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class S3BashLibTest extends Specification {

    // -- legacy

    def 'should get uploader script' () {

        given:
        def opts = Mock(AwsOptions)

        when:
        def script = S3BashLib.script(opts)
        then:
        1 * opts.getAwsCli() >> 'aws'
        1 * opts.getStorageClass() >> null
        1 * opts.getStorageEncryption() >> null

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
                    '''
                .stripIndent(true)
    }

    def 'should set storage class and encryption' () {

        given:
        def opts = Mock(AwsOptions)

        when:
        def script = S3BashLib.script(opts)
        then:
        opts.getStorageClass() >> 'S-CLAZZ'
        opts.getStorageEncryption() >> 'S-ENCRYPT'
        opts.getAwsCli() >> '/foo/bin/aws'
        opts.getMaxParallelTransfers() >> 33

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
                        IFS=$'\\n\'
                        local cmd=("$@")
                        local cpus=$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
                        local max=$(if (( cpus>33 )); then echo 33; else echo $cpus; fi)
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
                          /foo/bin/aws s3 cp --only-show-errors --sse S-ENCRYPT --storage-class S-CLAZZ - "$s3path"
                        elif [[ -d "$name" ]]; then
                          /foo/bin/aws s3 cp --only-show-errors --recursive --sse S-ENCRYPT --storage-class S-CLAZZ "$name" "$s3path/$name"
                        else
                          /foo/bin/aws s3 cp --only-show-errors --sse S-ENCRYPT --storage-class S-CLAZZ "$name" "$s3path/$name"
                        fi
                    }
                    
                    nxf_s3_download() {
                        local source=$1
                        local target=$2
                        local file_name=$(basename $1)
                        local is_dir=$(/foo/bin/aws s3 ls $source | grep -F "PRE ${file_name}/" -c)
                        if [[ $is_dir == 1 ]]; then
                            /foo/bin/aws s3 cp --only-show-errors --recursive "$source" "$target"
                        else 
                            /foo/bin/aws s3 cp --only-show-errors "$source" "$target"
                        fi
                    }
                    '''
                .stripIndent(true)

    }

    // -- new test

    def 'should create base script' () {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [:]
        }

        expect:
        S3BashLib.script()  == '''
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

    def 'should create base script with legacy retry mode' () {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [aws:[batch: [maxTransferAttempts: 100, retryMode: 'legacy']]]
        }

        expect:
        S3BashLib.script()  == '''
            # aws cli retry config
            export AWS_RETRY_MODE=legacy 
            export AWS_MAX_ATTEMPTS=100
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

    def 'should create base script with built-in retry mode' () {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [aws:[batch: [retryMode: 'built-in']]]
        }

        expect:
        S3BashLib.script()  == '''
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

    def 'should create base script with custom settings' () {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [aws:[batch: [cliPath: '/some/bin/aws', retryMode: 'legacy', maxTransferAttempts: 99]]]
        }

        expect:
        S3BashLib.script() == '''
            # aws cli retry config
            export AWS_RETRY_MODE=legacy 
            export AWS_MAX_ATTEMPTS=99
            # aws helper
            nxf_s3_upload() {
                local name=$1
                local s3path=$2
                if [[ "$name" == - ]]; then
                  /some/bin/aws s3 cp --only-show-errors --storage-class STANDARD - "$s3path"
                elif [[ -d "$name" ]]; then
                  /some/bin/aws s3 cp --only-show-errors --recursive --storage-class STANDARD "$name" "$s3path/$name"
                else
                  /some/bin/aws s3 cp --only-show-errors --storage-class STANDARD "$name" "$s3path/$name"
                fi
            }
            
            nxf_s3_download() {
                local source=$1
                local target=$2
                local file_name=$(basename $1)
                local is_dir=$(/some/bin/aws s3 ls $source | grep -F "PRE ${file_name}/" -c)
                if [[ $is_dir == 1 ]]; then
                    /some/bin/aws s3 cp --only-show-errors --recursive "$source" "$target"
                else 
                    /some/bin/aws s3 cp --only-show-errors "$source" "$target"
                fi
            }
            '''.stripIndent(true)
    }

    def 'should create base script with custom options' () {
        given:
        def opts = Mock(AwsOptions) {
            getMaxParallelTransfers() >> 5
            getMaxTransferAttempts() >> 10
            getDelayBetweenAttempts() >> Duration.of('20s')
        }

        expect:
        S3BashLib.script(opts)  == '''\
            # bash helper functions
            nxf_cp_retry() {
                local max_attempts=10
                local timeout=20
                local attempt=0
                local exitCode=0
                while (( $attempt < $max_attempts ))
                do
                  if "$@"
                    then
                      return 0
                  else
                    exitCode=$?
                  fi
                  if [[ $exitCode == 0 ]]
                  then
                    break
                  fi
                  nxf_sleep $timeout
                  attempt=$(( attempt + 1 ))
                  timeout=$(( timeout * 2 ))
                done
            }
            
            nxf_parallel() {
                IFS=$'\\n'
                local cmd=("$@")
                local cpus=$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
                local max=$(if (( cpus>5 )); then echo 5; else echo $cpus; fi)
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


    def 'should create base script with options' () {
        given:
        def opts = Mock(AwsOptions)

        expect:
        S3BashLib.script(opts)  == '''\
            # bash helper functions
            nxf_cp_retry() {
                local max_attempts=1
                local timeout=10
                local attempt=0
                local exitCode=0
                while (( $attempt < $max_attempts ))
                do
                  if "$@"
                    then
                      return 0
                  else
                    exitCode=$?
                  fi
                  if [[ $exitCode == 0 ]]
                  then
                    break
                  fi
                  nxf_sleep $timeout
                  attempt=$(( attempt + 1 ))
                  timeout=$(( timeout * 2 ))
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
    }

    def 'should create with storage encrypt' () {
        given:
        def sess1 = Mock(Session)  {
            getConfig() >> [aws: [ client: [ storageKmsKeyId: 'my-kms-key', storageEncryption: 'aws:kms']]]
        }
        and:
        def opts = new AwsOptions(sess1)

        expect:
        S3BashLib.script(opts)  == '''\
            # bash helper functions
            nxf_cp_retry() {
                local max_attempts=5
                local timeout=10
                local attempt=0
                local exitCode=0
                while (( $attempt < $max_attempts ))
                do
                  if "$@"
                    then
                      return 0
                  else
                    exitCode=$?
                  fi
                  if [[ $exitCode == 0 ]]
                  then
                    break
                  fi
                  nxf_sleep $timeout
                  attempt=$(( attempt + 1 ))
                  timeout=$(( timeout * 2 ))
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
            
            # aws cli retry config
            export AWS_RETRY_MODE=standard 
            export AWS_MAX_ATTEMPTS=5
            # aws helper
            nxf_s3_upload() {
                local name=$1
                local s3path=$2
                if [[ "$name" == - ]]; then
                  aws s3 cp --only-show-errors --sse aws:kms --sse-kms-key-id my-kms-key --storage-class STANDARD - "$s3path"
                elif [[ -d "$name" ]]; then
                  aws s3 cp --only-show-errors --recursive --sse aws:kms --sse-kms-key-id my-kms-key --storage-class STANDARD "$name" "$s3path/$name"
                else
                  aws s3 cp --only-show-errors --sse aws:kms --sse-kms-key-id my-kms-key --storage-class STANDARD "$name" "$s3path/$name"
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


    def 'should create with s3 acl' () {
        given:
        def sess1 = Mock(Session)  {
            getConfig() >> [aws: [ client: [ s3Acl: 'PublicRead']]]
        }
        and:
        def opts = new AwsOptions(sess1)

        expect:
        S3BashLib.script(opts)  == '''\
            # bash helper functions
            nxf_cp_retry() {
                local max_attempts=5
                local timeout=10
                local attempt=0
                local exitCode=0
                while (( $attempt < $max_attempts ))
                do
                  if "$@"
                    then
                      return 0
                  else
                    exitCode=$?
                  fi
                  if [[ $exitCode == 0 ]]
                  then
                    break
                  fi
                  nxf_sleep $timeout
                  attempt=$(( attempt + 1 ))
                  timeout=$(( timeout * 2 ))
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
            
            # aws cli retry config
            export AWS_RETRY_MODE=standard 
            export AWS_MAX_ATTEMPTS=5
            # aws helper
            nxf_s3_upload() {
                local name=$1
                local s3path=$2
                if [[ "$name" == - ]]; then
                  aws s3 cp --only-show-errors --acl public-read --storage-class STANDARD - "$s3path"
                elif [[ -d "$name" ]]; then
                  aws s3 cp --only-show-errors --recursive --acl public-read --storage-class STANDARD "$name" "$s3path/$name"
                else
                  aws s3 cp --only-show-errors --acl public-read --storage-class STANDARD "$name" "$s3path/$name"
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

    def 'should create s5cmd script' () {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [aws:[batch:[platformType: 'fargate', cliPath: 's5cmd']]]
        }

        expect:
        S3BashLib.script()  == '''
            # aws helper for s5cmd
            nxf_s3_upload() {
                local name=$1
                local s3path=$2
                if [[ "$name" == - ]]; then
                  local tmp=$(nxf_mktemp)
                  cp /dev/stdin $tmp/$name
                  s5cmd cp --storage-class STANDARD $tmp/$name "$s3path"
                elif [[ -d "$name" ]]; then
                  s5cmd cp --storage-class STANDARD "$name/" "$s3path/$name/"
                else
                  s5cmd cp --storage-class STANDARD "$name" "$s3path/$name"
                fi
            }
            
            nxf_s3_download() {
                local source=$1
                local target=$2
                local file_name=$(basename $1)
                local is_dir=$(s5cmd ls $source | grep -F "DIR  ${file_name}/" -c)
                if [[ $is_dir == 1 ]]; then
                    s5cmd cp "$source/*" "$target"
                else 
                    s5cmd cp "$source" "$target"
                fi
            }
            '''.stripIndent(true)
    }

    def 'should create s5cmd script with acl' () {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [aws:[batch:[platformType: 'fargate', cliPath: 's5cmd'], client:[ s3Acl: 'PublicRead']]]
        }

        expect:
        S3BashLib.script()  == '''
            # aws helper for s5cmd
            nxf_s3_upload() {
                local name=$1
                local s3path=$2
                if [[ "$name" == - ]]; then
                  local tmp=$(nxf_mktemp)
                  cp /dev/stdin $tmp/$name
                  s5cmd cp --acl public-read --storage-class STANDARD $tmp/$name "$s3path"
                elif [[ -d "$name" ]]; then
                  s5cmd cp --acl public-read --storage-class STANDARD "$name/" "$s3path/$name/"
                else
                  s5cmd cp --acl public-read --storage-class STANDARD "$name" "$s3path/$name"
                fi
            }
            
            nxf_s3_download() {
                local source=$1
                local target=$2
                local file_name=$(basename $1)
                local is_dir=$(s5cmd ls $source | grep -F "DIR  ${file_name}/" -c)
                if [[ $is_dir == 1 ]]; then
                    s5cmd cp "$source/*" "$target"
                else 
                    s5cmd cp "$source" "$target"
                fi
            }
            '''.stripIndent(true)
    }
    

}
