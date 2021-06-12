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

    def 'should create base script' () {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [:]
        }

        expect:
        S3BashLib.script()  == '''
            # aws helper
            nxf_s3_upload() {
                local name=$1
                local s3path=$2
                if [[ -d "$name" ]]; then
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
            '''.stripIndent()
    }

    def 'should create base script with custom settings' () {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [aws:[batch: [cliPath: '/some/bin/aws', retryMode: 'legacy', maxTransferAttempts: 99]]]
        }

        expect:
        S3BashLib.script()  == '''
            # aws cli retry config
            export AWS_RETRY_MODE=legacy 
            export AWS_MAX_ATTEMPTS=99
            # aws helper
            nxf_s3_upload() {
                local name=$1
                local s3path=$2
                if [[ -d "$name" ]]; then
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
            '''.stripIndent()
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
                      [[ -e /proc/$x ]] && copy+=($x)
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
                if [[ -d "$name" ]]; then
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
            '''.stripIndent()
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
                      [[ -e /proc/$x ]] && copy+=($x)
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
                if [[ -d "$name" ]]; then
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
            '''.stripIndent()
    }

}
