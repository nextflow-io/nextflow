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
import java.nio.file.Path
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

            nxf_az_download() {
                local source=$1
                local target=$2
                local basedir=$(dirname "$2")
                local ret
                mkdir -p "$basedir"

                azcopy cp "$source?$AZ_SAS" "$target" >/dev/null 2>&1 || {
                    ret=$?
                    ## killed by a signal (e.g. out of memory), retrying as a directory would hide the cause
                    if (( ret > 128 )); then
                        >&2 echo "Unable to download path: $source (azcopy exit status $ret)"
                        return $ret
                    fi
                    ## if fails check if it was trying to download a directory
                    mkdir -p "$target"
                    azcopy cp "$source/*?$AZ_SAS" "$target" --recursive >/dev/null || {
                        ret=$?
                        rm -rf "$target"
                        >&2 echo "Unable to download path: $source (azcopy exit status $ret)"
                        return $ret
                    }
                }
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

            nxf_az_download() {
                local source=$1
                local target=$2
                local basedir=$(dirname "$2")
                local ret
                mkdir -p "$basedir"

                azcopy cp "$source?$AZ_SAS" "$target" >/dev/null 2>&1 || {
                    ret=$?
                    ## killed by a signal (e.g. out of memory), retrying as a directory would hide the cause
                    if (( ret > 128 )); then
                        >&2 echo "Unable to download path: $source (azcopy exit status $ret)"
                        return $ret
                    fi
                    ## if fails check if it was trying to download a directory
                    mkdir -p "$target"
                    azcopy cp "$source/*?$AZ_SAS" "$target" --recursive >/dev/null || {
                        ret=$?
                        rm -rf "$target"
                        >&2 echo "Unable to download path: $source (azcopy exit status $ret)"
                        return $ret
                    }
                }
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

            nxf_az_download() {
                local source=$1
                local target=$2
                local basedir=$(dirname "$2")
                local ret
                mkdir -p "$basedir"

                azcopy cp "$source?$AZ_SAS" "$target" >/dev/null 2>&1 || {
                    ret=$?
                    ## killed by a signal (e.g. out of memory), retrying as a directory would hide the cause
                    if (( ret > 128 )); then
                        >&2 echo "Unable to download path: $source (azcopy exit status $ret)"
                        return $ret
                    fi
                    ## if fails check if it was trying to download a directory
                    mkdir -p "$target"
                    azcopy cp "$source/*?$AZ_SAS" "$target" --recursive >/dev/null || {
                        ret=$?
                        rm -rf "$target"
                        >&2 echo "Unable to download path: $source (azcopy exit status $ret)"
                        return $ret
                    }
                }
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
        // each attempt tries the file copy and then the directory copy
        workDir.resolve('calls').readLines().size() == 6

        cleanup:
        workDir?.toFile()?.deleteDir()
    }

    private static final String FAKE_AZCOPY = '''\
        #!/bin/bash
        echo "$2" >> calls
        src=${2%%\\?*}
        dst=$3
        case $src in
          */file.txt)    echo data > "$dst" ;;
          */oom.txt)     echo partial > "$dst"; kill -KILL $$ ;;
          */dir/\\*)      echo a > "$dst/a.txt" ;;
          */oomdir/\\*)   echo a > "$dst/a.txt"; kill -KILL $$ ;;
          *)             echo "failed to perform copy command due to error: cannot use directory as source without --recursive or a trailing wildcard (/*)"; exit 1 ;;
        esac
        '''.stripIndent()

    private Path createFakeAzcopy(Path workDir) {
        def bin = Files.createDirectory(workDir.resolve('bin'))
        def azcopy = bin.resolve('azcopy')
        azcopy.text = FAKE_AZCOPY
        azcopy.toFile().setExecutable(true)
        return bin
    }

    @Unroll
    def 'should download #SOURCE and return #STATUS' () {
        given:
        def workDir = Files.createTempDirectory('az-download')
        def bin = createFakeAzcopy(workDir)
        def script = workDir.resolve('test.sh')
        // run in the background, as nxf_parallel does, so that errexit applies inside the function
        script.text = 'set -eu\n' + AzBashLib.script() + """
            export AZ_SAS=token
            export PATH=${bin}:\$PATH
            nxf_az_download 'https://acct.blob.core.windows.net/c/${SOURCE}' target 2> stderr &
            wait \$! && echo 0 > status || echo \$? > status
            """.stripIndent()

        when:
        def process = ['bash', script.toString()].execute(null as List, workDir.toFile())
        process.waitFor(30, TimeUnit.SECONDS)
        def target = workDir.resolve('target')
        def stderr = workDir.resolve('stderr').text

        then:
        process.exitValue() == 0
        workDir.resolve('status').text.trim() as int == STATUS
        (Files.isDirectory(target) ? 'dir' : Files.exists(target) ? 'file' : 'none') == RESULT
        workDir.resolve('calls').readLines().size() == CALLS
        ERROR ? stderr.endsWith(ERROR) : stderr.isEmpty()
        !stderr.contains('File exists')

        cleanup:
        workDir?.toFile()?.deleteDir()

        where:
        SOURCE        | STATUS | RESULT | CALLS | ERROR
        'file.txt'    | 0      | 'file' | 1     | ''
        'dir'         | 0      | 'dir'  | 2     | ''
        // the real azcopy exits 0 on a recursive copy of a missing path, leaving an empty directory
        'missing'     | 1      | 'none' | 2     | 'Unable to download path: https://acct.blob.core.windows.net/c/missing (azcopy exit status 1)\n'
        'oom.txt'     | 137    | 'file' | 1     | 'Unable to download path: https://acct.blob.core.windows.net/c/oom.txt (azcopy exit status 137)\n'
        'oomdir'      | 137    | 'none' | 2     | 'Unable to download path: https://acct.blob.core.windows.net/c/oomdir (azcopy exit status 137)\n'
    }

}
