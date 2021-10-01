package nextflow.cloud.azure.file

import nextflow.util.Duration
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzBashLibTest extends Specification{

    def 'should return base script' () {
        expect:
        AzBashLib.script() == '''
                nxf_az_upload() {
                    local name=$1
                    local target=${2%/} ## remove ending slash
                
                    if [[ -d $name ]]; then
                      azcopy cp "$name" "$target?$AZ_SAS" --recursive
                    else
                      azcopy cp "$name" "$target/$name?$AZ_SAS"
                    fi
                }
                
                nxf_az_download() {
                    local source=$1
                    local target=$2
                    local basedir=$(dirname $2)
                    local ret
                    mkdir -p "$basedir"
                
                    ret=$(azcopy cp "$source?$AZ_SAS" "$target" 2>&1) || {
                        ## if fails check if it was trying to download a directory
                        mkdir -p $target
                        azcopy cp "$source/*?$AZ_SAS" "$target" --recursive >/dev/null || {
                            rm -rf $target
                            >&2 echo "Unable to download path: $source"
                            exit 1
                        }
                    }
                }
                '''.stripIndent()
    }

    def 'should return script with config' () {
        expect:
        AzBashLib.script(10, 20, Duration.of('30s')) == '''\
            # bash helper functions
            nxf_cp_retry() {
                local max_attempts=20
                local timeout=30
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
                local max=$(if (( cpus>10 )); then echo 10; else echo $cpus; fi)
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
            
            nxf_az_upload() {
                local name=$1
                local target=${2%/} ## remove ending slash
            
                if [[ -d $name ]]; then
                  azcopy cp "$name" "$target?$AZ_SAS" --recursive
                else
                  azcopy cp "$name" "$target/$name?$AZ_SAS"
                fi
            }
            
            nxf_az_download() {
                local source=$1
                local target=$2
                local basedir=$(dirname $2)
                local ret
                mkdir -p "$basedir"
            
                ret=$(azcopy cp "$source?$AZ_SAS" "$target" 2>&1) || {
                    ## if fails check if it was trying to download a directory
                    mkdir -p $target
                    azcopy cp "$source/*?$AZ_SAS" "$target" --recursive >/dev/null || {
                        rm -rf $target
                        >&2 echo "Unable to download path: $source"
                        exit 1
                    }
                }
            }
            '''.stripIndent()
    }
}
