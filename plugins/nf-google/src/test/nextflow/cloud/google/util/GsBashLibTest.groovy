package nextflow.cloud.google.util


import nextflow.Session
import nextflow.cloud.google.lifesciences.GoogleLifeSciencesConfig
import nextflow.util.Duration
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GsBashLibTest extends Specification {

    def 'should make gs env' () {
        expect:
        new GsBashLib().makeEnv() == '''
                # google storage helper
                gs_opts=('-q' '-m' '-o' 'GSUtil:parallel_thread_count=1' '-o' 'GSUtil:sliced_object_download_max_components=8')
                '''.stripIndent()

        and:
        new GsBashLib()
                .withParallelThreadCount(2)
                .withDownloadMaxComponents(4)
                .makeEnv() == '''
                # google storage helper
                gs_opts=('-q' '-m' '-o' 'GSUtil:parallel_thread_count=2' '-o' 'GSUtil:sliced_object_download_max_components=4')
                '''.stripIndent()
    }

    def 'should make lib'() {
        expect:
        new GsBashLib().makeLib() == '''
            nxf_gs_download() {
                local source=$1
                local target=$2
                local project=$3
                local basedir=$(dirname $2)
                local ret
                local opts=("${gs_opts[@]}")
                if [[ $project ]]; then
                  opts+=('-u' "$project")
                fi
            
                ## download assuming it's a file download
                mkdir -p $basedir
                ret=$(gsutil ${opts[@]} cp "$source" "$target" 2>&1) || {
                    ## if fails check if it was trying to download a directory
                    mkdir $target
                    gsutil ${opts[@]} cp -R "$source/*" "$target" || {
                      rm -rf $target
                      >&2 echo "Unable to download path: $source"
                      exit 1
                    }
                }
            }
            
            nxf_gs_upload() {
                local name=$1
                local target=$2
                gsutil ${gs_opts[@]} cp -R "$name" "$target/$name"
            }
            '''.stripIndent()
    }

    def 'should make from session' () {
        when:
        def session = Mock(Session) {
            getConfig() >> [google: [storage:[parallelThreadCount: 3, downloadMaxComponents: 9]]]
        }
        and:
        def result = GsBashLib.fromSession(session)
        then:
        result == '''
            # google storage helper
            gs_opts=('-q' '-m' '-o' 'GSUtil:parallel_thread_count=3' '-o' 'GSUtil:sliced_object_download_max_components=9')
            
            nxf_gs_download() {
                local source=$1
                local target=$2
                local project=$3
                local basedir=$(dirname $2)
                local ret
                local opts=("${gs_opts[@]}")
                if [[ $project ]]; then
                  opts+=('-u' "$project")
                fi
            
                ## download assuming it's a file download
                mkdir -p $basedir
                ret=$(gsutil ${opts[@]} cp "$source" "$target" 2>&1) || {
                    ## if fails check if it was trying to download a directory
                    mkdir $target
                    gsutil ${opts[@]} cp -R "$source/*" "$target" || {
                      rm -rf $target
                      >&2 echo "Unable to download path: $source"
                      exit 1
                    }
                }
            }
            
            nxf_gs_upload() {
                local name=$1
                local target=$2
                gsutil ${gs_opts[@]} cp -R "$name" "$target/$name"
            }
            '''.stripIndent()
    }

    def 'should make from config' () {
        when:
        def config = Mock(GoogleLifeSciencesConfig)

        and:
        def result = GsBashLib.fromConfig(config)
        then:
        result == '''\
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
            
            # google storage helper
            gs_opts=('-q' '-m' '-o' 'GSUtil:parallel_thread_count=1' '-o' 'GSUtil:sliced_object_download_max_components=8')
            
            nxf_gs_download() {
                local source=$1
                local target=$2
                local project=$3
                local basedir=$(dirname $2)
                local ret
                local opts=("${gs_opts[@]}")
                if [[ $project ]]; then
                  opts+=('-u' "$project")
                fi
            
                ## download assuming it's a file download
                mkdir -p $basedir
                ret=$(gsutil ${opts[@]} cp "$source" "$target" 2>&1) || {
                    ## if fails check if it was trying to download a directory
                    mkdir $target
                    gsutil ${opts[@]} cp -R "$source/*" "$target" || {
                      rm -rf $target
                      >&2 echo "Unable to download path: $source"
                      exit 1
                    }
                }
            }
            
            nxf_gs_upload() {
                local name=$1
                local target=$2
                gsutil ${gs_opts[@]} cp -R "$name" "$target/$name"
            }
            '''.stripIndent()
    }


    def 'should make from config and custom settings' () {
        when:
        def config = Mock(GoogleLifeSciencesConfig) {
            getMaxParallelTransfers() >> 5
            getMaxTransferAttempts() >> 10
            getDelayBetweenAttempts() >> Duration.of('20s')
            getParallelThreadCount() >> 20
            getDownloadMaxComponents() >> 30
        }

        and:
        def result = GsBashLib.fromConfig(config)
        then:
        result == '''\
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
            
            # google storage helper
            gs_opts=('-q' '-m' '-o' 'GSUtil:parallel_thread_count=20' '-o' 'GSUtil:sliced_object_download_max_components=30')
            
            nxf_gs_download() {
                local source=$1
                local target=$2
                local project=$3
                local basedir=$(dirname $2)
                local ret
                local opts=("${gs_opts[@]}")
                if [[ $project ]]; then
                  opts+=('-u' "$project")
                fi
            
                ## download assuming it's a file download
                mkdir -p $basedir
                ret=$(gsutil ${opts[@]} cp "$source" "$target" 2>&1) || {
                    ## if fails check if it was trying to download a directory
                    mkdir $target
                    gsutil ${opts[@]} cp -R "$source/*" "$target" || {
                      rm -rf $target
                      >&2 echo "Unable to download path: $source"
                      exit 1
                    }
                }
            }
            
            nxf_gs_upload() {
                local name=$1
                local target=$2
                gsutil ${gs_opts[@]} cp -R "$name" "$target/$name"
            }
            '''.stripIndent()
    }
}
