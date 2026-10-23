/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cloud.azure.file


import java.nio.file.Files
import java.util.concurrent.TimeUnit

import nextflow.cloud.azure.config.AzCopyOpts
import nextflow.util.Duration
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzBashLibTest extends Specification {

    def 'should return base script'() {
        expect:
        AzBashLib.script() == '''
            # custom env variables used for azcopy opts
            export AZCOPY_BLOCK_SIZE_MB=4
            export AZCOPY_BLOCK_BLOB_TIER=None

            nxf_az_upload() {
                local name=$1
                local target=${2%/} ## remove ending slash
                local base_name="$(basename "$name")"
                local dir_name="$(dirname "$name")"

                if [[ -d $name ]]; then
                  if [[ "$base_name" == "$name" ]]; then
                    azcopy cp "$name" "$target?$AZ_SAS" --recursive --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
                  else
                    azcopy cp "$name" "$target/$dir_name?$AZ_SAS" --recursive --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
                  fi
                else
                  azcopy cp "$name" "$target/$name?$AZ_SAS" --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
                fi
            }

            ## prints the last progress update and the last messages of the azcopy output
            nxf_az_log() {
                printf '%s\\n' "$1" | tr '\\r' '\\n' | sed 's/[[:space:]]*$//' | grep -v '^$' | awk '/^[0-9.]+ %, /{p=$0; next} {l[n++]=$0} END {if(p!="") print p; for(i=(n>10?n-10:0);i<n;i++) print l[i]}' >&2 || true
            }

            nxf_az_download() {
                local source=$1
                local target=$2
                local basedir=$(dirname $2)
                local out
                local ret=0
                mkdir -p "$basedir"

                out=$(azcopy cp "$source?$AZ_SAS" "$target" 2>&1) || ret=$?
                if [[ $ret != 0 ]]; then
                    >&2 echo "Unable to download file: $source"
                    nxf_az_log "$out"
                    return $ret
                fi
            }

            nxf_az_download_dir() {
                local source=$1
                local target=$2
                local out
                local ret=0
                mkdir -p "$target"

                out=$(azcopy cp "$source/*?$AZ_SAS" "$target" --recursive 2>&1) || ret=$?
                if [[ $ret != 0 ]]; then
                    rm -rf "$target"
                    >&2 echo "Unable to download directory: $source"
                    nxf_az_log "$out"
                    return $ret
                fi
            }
            '''.stripIndent()
    }

    def 'should return script with config, with default azcopy opts'() {

        expect:
        AzBashLib.script(new AzCopyOpts(), 10, 20, Duration.of('30s')) == '''\
            # bash helper functions
            nxf_cp_retry() {
                local max_attempts=20
                local timeout=30
                local max_timeout=300
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
                  attempt=$(( attempt + 1 ))
                  if (( $attempt >= $max_attempts || $exitCode > 128 ))
                  then
                    break
                  fi
                  nxf_sleep $timeout
                  timeout=$(( timeout * 2 > max_timeout ? max_timeout : timeout * 2 ))
                done
                return $exitCode
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

            # custom env variables used for azcopy opts
            export AZCOPY_BLOCK_SIZE_MB=4
            export AZCOPY_BLOCK_BLOB_TIER=None

            nxf_az_upload() {
                local name=$1
                local target=${2%/} ## remove ending slash
                local base_name="$(basename "$name")"
                local dir_name="$(dirname "$name")"

                if [[ -d $name ]]; then
                  if [[ "$base_name" == "$name" ]]; then
                    azcopy cp "$name" "$target?$AZ_SAS" --recursive --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
                  else
                    azcopy cp "$name" "$target/$dir_name?$AZ_SAS" --recursive --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
                  fi
                else
                  azcopy cp "$name" "$target/$name?$AZ_SAS" --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
                fi
            }

            ## prints the last progress update and the last messages of the azcopy output
            nxf_az_log() {
                printf '%s\\n' "$1" | tr '\\r' '\\n' | sed 's/[[:space:]]*$//' | grep -v '^$' | awk '/^[0-9.]+ %, /{p=$0; next} {l[n++]=$0} END {if(p!="") print p; for(i=(n>10?n-10:0);i<n;i++) print l[i]}' >&2 || true
            }

            nxf_az_download() {
                local source=$1
                local target=$2
                local basedir=$(dirname $2)
                local out
                local ret=0
                mkdir -p "$basedir"

                out=$(azcopy cp "$source?$AZ_SAS" "$target" 2>&1) || ret=$?
                if [[ $ret != 0 ]]; then
                    >&2 echo "Unable to download file: $source"
                    nxf_az_log "$out"
                    return $ret
                fi
            }

            nxf_az_download_dir() {
                local source=$1
                local target=$2
                local out
                local ret=0
                mkdir -p "$target"

                out=$(azcopy cp "$source/*?$AZ_SAS" "$target" --recursive 2>&1) || ret=$?
                if [[ $ret != 0 ]]; then
                    rm -rf "$target"
                    >&2 echo "Unable to download directory: $source"
                    nxf_az_log "$out"
                    return $ret
                fi
            }
            '''.stripIndent()
    }

    def 'should return script with config, with custom azcopy opts'() {

        expect:
        AzBashLib.script(new AzCopyOpts([blobTier: "Hot", blockSize: "10"]), 10, 20, Duration.of('30s')) == '''\
            # bash helper functions
            nxf_cp_retry() {
                local max_attempts=20
                local timeout=30
                local max_timeout=300
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
                  attempt=$(( attempt + 1 ))
                  if (( $attempt >= $max_attempts || $exitCode > 128 ))
                  then
                    break
                  fi
                  nxf_sleep $timeout
                  timeout=$(( timeout * 2 > max_timeout ? max_timeout : timeout * 2 ))
                done
                return $exitCode
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

            # custom env variables used for azcopy opts
            export AZCOPY_BLOCK_SIZE_MB=10
            export AZCOPY_BLOCK_BLOB_TIER=Hot

            nxf_az_upload() {
                local name=$1
                local target=${2%/} ## remove ending slash
                local base_name="$(basename "$name")"
                local dir_name="$(dirname "$name")"

                if [[ -d $name ]]; then
                  if [[ "$base_name" == "$name" ]]; then
                    azcopy cp "$name" "$target?$AZ_SAS" --recursive --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
                  else
                    azcopy cp "$name" "$target/$dir_name?$AZ_SAS" --recursive --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
                  fi
                else
                  azcopy cp "$name" "$target/$name?$AZ_SAS" --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
                fi
            }

            ## prints the last progress update and the last messages of the azcopy output
            nxf_az_log() {
                printf '%s\\n' "$1" | tr '\\r' '\\n' | sed 's/[[:space:]]*$//' | grep -v '^$' | awk '/^[0-9.]+ %, /{p=$0; next} {l[n++]=$0} END {if(p!="") print p; for(i=(n>10?n-10:0);i<n;i++) print l[i]}' >&2 || true
            }

            nxf_az_download() {
                local source=$1
                local target=$2
                local basedir=$(dirname $2)
                local out
                local ret=0
                mkdir -p "$basedir"

                out=$(azcopy cp "$source?$AZ_SAS" "$target" 2>&1) || ret=$?
                if [[ $ret != 0 ]]; then
                    >&2 echo "Unable to download file: $source"
                    nxf_az_log "$out"
                    return $ret
                fi
            }

            nxf_az_download_dir() {
                local source=$1
                local target=$2
                local out
                local ret=0
                mkdir -p "$target"

                out=$(azcopy cp "$source/*?$AZ_SAS" "$target" --recursive 2>&1) || ret=$?
                if [[ $ret != 0 ]]; then
                    rm -rf "$target"
                    >&2 echo "Unable to download directory: $source"
                    nxf_az_log "$out"
                    return $ret
                fi
            }
            '''.stripIndent()
    }

    def 'should retry failed downloads' () {
        given:
        def workDir = Files.createTempDirectory('az-retry')
        def bin = Files.createDirectory(workDir.resolve('bin'))
        def azcopy = bin.resolve('azcopy')
        azcopy.text = '#!/bin/bash\necho x >> calls\nexit 1\n'
        azcopy.toFile().setExecutable(true)
        def script = workDir.resolve('test.sh')
        script.text = AzBashLib.script(new AzCopyOpts(), 4, 3, Duration.of('1 sec')) + """
            nxf_sleep() { :; }
            export PATH=${bin}:\$PATH
            nxf_cp_retry nxf_az_download 'https://acct.blob.core.windows.net/c/file.txt' file.txt 2>/dev/null && echo 0 > status || echo \$? > status
            """.stripIndent()

        when:
        def process = ['bash', script.toString()].execute(null as List, workDir.toFile())
        process.waitFor(30, TimeUnit.SECONDS)

        then:
        process.exitValue() == 0
        workDir.resolve('status').text.trim() == '1'
        workDir.resolve('calls').readLines().size() == 3

        cleanup:
        workDir?.toFile()?.deleteDir()
    }

    @Unroll
    def 'should #DOWNLOAD #SOURCE and return #STATUS' () {
        given:
        def workDir = Files.createTempDirectory('az-download')
        def bin = Files.createDirectory(workDir.resolve('bin'))
        def azcopy = bin.resolve('azcopy')
        // mimics azcopy 10.8: a non-recursive copy of a missing blob fails as if it were a
        // directory, and a recursive copy of a missing or empty prefix transfers nothing
        azcopy.text = """\
            #!/bin/bash
            src=\${2%%\\?*}
            dst=\$3
            case \$src in
              */file.txt)       echo data > \$dst ;;
              */oom.txt)        echo partial > \$dst; for i in \$(seq 1 500); do printf '%s %%, 0 Done, 0 Failed, 1 Pending, 0 Skipped, 1 Total, 2-sec Throughput (Mb/s): 1\\r' \$i; done; exit 137 ;;
              */denied.txt)     echo "X-Ms-Error-Code: [AuthenticationFailed]"; exit 1 ;;
              */dir/\\*)        echo a > \$dst/a.txt ;;
              */marker/\\*)     touch \$dst/.azure_blob_dir ;;
              */oomdir/\\*)     echo a > \$dst/a.txt; exit 137 ;;
              */denied/\\*)     echo "X-Ms-Error-Code: [AuthenticationFailed]"; exit 1 ;;
              */\\*)            echo "Total Number of Transfers: 0" ;;
              *)                echo "failed to perform copy command due to error: cannot use directory as source without --recursive or a trailing wildcard (/*)"; exit 1 ;;
            esac
            """.stripIndent()
        azcopy.toFile().setExecutable(true)
        def script = workDir.resolve('test.sh')
        script.text = 'set -eu\n' + AzBashLib.script() + """
            export AZ_SAS=token
            export PATH=${bin}:\$PATH
            ${DOWNLOAD} 'https://acct.blob.core.windows.net/c/${SOURCE}' target 2> stderr && echo 0 > status || echo \$? > status
            """.stripIndent()

        when:
        def process = ['bash', script.toString()].execute(null as List, workDir.toFile())
        process.waitFor(30, TimeUnit.SECONDS)
        def target = workDir.resolve('target')
        def stderr = workDir.resolve('stderr').text

        then:
        process.exitValue() == 0
        workDir.resolve('status').text.trim() as int == STATUS
        (Files.isDirectory(target) ? "dir:${target.toFile().list().sort().join(',')}" : Files.exists(target) ? 'file' : 'none') == RESULT
        stderr.contains(ERROR)
        !stderr.contains('mkdir')
        stderr.size() < 1000

        cleanup:
        workDir?.toFile()?.deleteDir()

        where:
        DOWNLOAD                | SOURCE        | STATUS | RESULT                | ERROR
        'nxf_az_download'       | 'file.txt'    | 0      | 'file'                | ''
        'nxf_az_download'       | 'missing.txt' | 1      | 'none'                | 'Unable to download file'
        'nxf_az_download'       | 'denied.txt'  | 1      | 'none'                | 'AuthenticationFailed'
        'nxf_az_download'       | 'oom.txt'     | 137    | 'file'                | 'Unable to download file'
        'nxf_az_download_dir'   | 'dir'         | 0      | 'dir:a.txt'           | ''
        'nxf_az_download_dir'   | 'marker'      | 0      | 'dir:.azure_blob_dir' | ''
        'nxf_az_download_dir'   | 'empty'       | 0      | 'dir:'                | ''
        'nxf_az_download_dir'   | 'denied'      | 1      | 'none'                | 'Unable to download directory'
        'nxf_az_download_dir'   | 'oomdir'      | 137    | 'none'                | 'Unable to download directory'
    }

}
