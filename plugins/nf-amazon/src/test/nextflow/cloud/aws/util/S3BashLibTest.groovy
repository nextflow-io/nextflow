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

    def 'should create base script with default retry mode' () {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [:] }

        expect:
        S3BashLib.script()  == '''
            # aws cli retry config
            export AWS_RETRY_MODE=standard 
            export AWS_MAX_ATTEMPTS=5
            # aws helper
            nxf_s3_upload() {
                local name=$1
                local s3path=$2
                shift 2
                # Collect remaining args in an array to preserve quoting & handle empty safely
                local opts=( "$@" )
                if [[ "$name" == - ]]; then
                  aws s3 cp --only-show-errors "${opts[@]}" - "$s3path"
                elif [[ -d "$name" ]]; then
                  aws s3 cp --only-show-errors --recursive "${opts[@]}" "$name" "$s3path/$name"
                else
                  aws s3 cp --only-show-errors "${opts[@]}" "$name" "$s3path/$name"
                fi
            }
            
            nxf_s3_download() {
                local source=$1
                local target=$2
                local file_name=$(basename $1)
                shift 2
                # Collect remaining args in an array to preserve quoting & handle empty safely
                local opts=( "$@" )
                local is_dir=$(aws s3 ls "${opts[@]}" $source | grep -F "PRE ${file_name}/" -c)
                if [[ $is_dir == 1 ]]; then
                    aws s3 cp --only-show-errors --recursive "${opts[@]}" "$source" "$target"
                else 
                    aws s3 cp --only-show-errors "${opts[@]}" "$source" "$target"
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
                shift 2
                # Collect remaining args in an array to preserve quoting & handle empty safely
                local opts=( "$@" )
                if [[ "$name" == - ]]; then
                  aws s3 cp --only-show-errors "${opts[@]}" - "$s3path"
                elif [[ -d "$name" ]]; then
                  aws s3 cp --only-show-errors --recursive "${opts[@]}" "$name" "$s3path/$name"
                else
                  aws s3 cp --only-show-errors "${opts[@]}" "$name" "$s3path/$name"
                fi
            }
            
            nxf_s3_download() {
                local source=$1
                local target=$2
                local file_name=$(basename $1)
                shift 2
                # Collect remaining args in an array to preserve quoting & handle empty safely
                local opts=( "$@" )
                local is_dir=$(aws s3 ls "${opts[@]}" $source | grep -F "PRE ${file_name}/" -c)
                if [[ $is_dir == 1 ]]; then
                    aws s3 cp --only-show-errors --recursive "${opts[@]}" "$source" "$target"
                else 
                    aws s3 cp --only-show-errors "${opts[@]}" "$source" "$target"
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
                shift 2
                # Collect remaining args in an array to preserve quoting & handle empty safely
                local opts=( "$@" )
                if [[ "$name" == - ]]; then
                  aws s3 cp --only-show-errors "${opts[@]}" - "$s3path"
                elif [[ -d "$name" ]]; then
                  aws s3 cp --only-show-errors --recursive "${opts[@]}" "$name" "$s3path/$name"
                else
                  aws s3 cp --only-show-errors "${opts[@]}" "$name" "$s3path/$name"
                fi
            }
            
            nxf_s3_download() {
                local source=$1
                local target=$2
                local file_name=$(basename $1)
                shift 2
                # Collect remaining args in an array to preserve quoting & handle empty safely
                local opts=( "$@" )
                local is_dir=$(aws s3 ls "${opts[@]}" $source | grep -F "PRE ${file_name}/" -c)
                if [[ $is_dir == 1 ]]; then
                    aws s3 cp --only-show-errors --recursive "${opts[@]}" "$source" "$target"
                else 
                    aws s3 cp --only-show-errors "${opts[@]}" "$source" "$target"
                fi
            }
            '''.stripIndent(true)
    }

    def 'should create base script with custom cli path' () {
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
                shift 2
                # Collect remaining args in an array to preserve quoting & handle empty safely
                local opts=( "$@" )
                if [[ "$name" == - ]]; then
                  /some/bin/aws s3 cp --only-show-errors "${opts[@]}" - "$s3path"
                elif [[ -d "$name" ]]; then
                  /some/bin/aws s3 cp --only-show-errors --recursive "${opts[@]}" "$name" "$s3path/$name"
                else
                  /some/bin/aws s3 cp --only-show-errors "${opts[@]}" "$name" "$s3path/$name"
                fi
            }
            
            nxf_s3_download() {
                local source=$1
                local target=$2
                local file_name=$(basename $1)
                shift 2
                # Collect remaining args in an array to preserve quoting & handle empty safely
                local opts=( "$@" )
                local is_dir=$(/some/bin/aws s3 ls "${opts[@]}" $source | grep -F "PRE ${file_name}/" -c)
                if [[ $is_dir == 1 ]]; then
                    /some/bin/aws s3 cp --only-show-errors --recursive "${opts[@]}" "$source" "$target"
                else 
                    /some/bin/aws s3 cp --only-show-errors "${opts[@]}" "$source" "$target"
                fi
            }
            '''.stripIndent(true)
    }

    def 'should create script with custom options including core functions' () {
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
                shift 2
                # Collect remaining args in an array to preserve quoting & handle empty safely
                local opts=( "$@" )
                if [[ "$name" == - ]]; then
                  aws s3 cp --only-show-errors "${opts[@]}" - "$s3path"
                elif [[ -d "$name" ]]; then
                  aws s3 cp --only-show-errors --recursive "${opts[@]}" "$name" "$s3path/$name"
                else
                  aws s3 cp --only-show-errors "${opts[@]}" "$name" "$s3path/$name"
                fi
            }
            
            nxf_s3_download() {
                local source=$1
                local target=$2
                local file_name=$(basename $1)
                shift 2
                # Collect remaining args in an array to preserve quoting & handle empty safely
                local opts=( "$@" )
                local is_dir=$(aws s3 ls "${opts[@]}" $source | grep -F "PRE ${file_name}/" -c)
                if [[ $is_dir == 1 ]]; then
                    aws s3 cp --only-show-errors --recursive "${opts[@]}" "$source" "$target"
                else 
                    aws s3 cp --only-show-errors "${opts[@]}" "$source" "$target"
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
                shift 2
                # Collect remaining args in an array to preserve quoting & handle empty safely
                local opts=( "$@" )
                if [[ "$name" == - ]]; then
                  aws s3 cp --only-show-errors "${opts[@]}" - "$s3path"
                elif [[ -d "$name" ]]; then
                  aws s3 cp --only-show-errors --recursive "${opts[@]}" "$name" "$s3path/$name"
                else
                  aws s3 cp --only-show-errors "${opts[@]}" "$name" "$s3path/$name"
                fi
            }
            
            nxf_s3_download() {
                local source=$1
                local target=$2
                local file_name=$(basename $1)
                shift 2
                # Collect remaining args in an array to preserve quoting & handle empty safely
                local opts=( "$@" )
                local is_dir=$(aws s3 ls "${opts[@]}" $source | grep -F "PRE ${file_name}/" -c)
                if [[ $is_dir == 1 ]]; then
                    aws s3 cp --only-show-errors --recursive "${opts[@]}" "$source" "$target"
                else 
                    aws s3 cp --only-show-errors "${opts[@]}" "$source" "$target"
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
                shift 2
                # Collect remaining args in an array to preserve quoting & handle empty safely
                local opts=( "$@" )
                if [[ "$name" == - ]]; then
                  local tmp=$(nxf_mktemp)
                  cp /dev/stdin $tmp/$name
                  s5cmd cp "${opts[@]}" $tmp/$name "$s3path"
                elif [[ -d "$name" ]]; then
                  s5cmd cp "${opts[@]}" "$name/" "$s3path/$name/"
                else
                  s5cmd cp "${opts[@]}" "$name" "$s3path/$name"
                fi
            }
            
            nxf_s3_download() {
                local source=$1
                local target=$2
                local file_name=$(basename $1)
                shift 2
                # Collect remaining args in an array to preserve quoting & handle empty safely
                local opts=( "$@" )
                local is_dir=$(s5cmd ls ${opts[@]} $source | grep -F "DIR  ${file_name}/" -c)
                if [[ $is_dir == 1 ]]; then
                    s5cmd cp ${opts[@]} "$source/*" "$target"
                else 
                    s5cmd cp ${opts[@]} "$source" "$target"
                fi
            }
            '''.stripIndent(true)
    }

}