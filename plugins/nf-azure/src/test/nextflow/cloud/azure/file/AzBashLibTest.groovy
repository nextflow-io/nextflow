package nextflow.cloud.azure.file


import nextflow.cloud.azure.config.AzCopyOpts
import nextflow.util.Duration
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzBashLibTest extends Specification {

    /**
     * Helper method to compare scripts while being tolerant of whitespace differences
     */
    boolean compareScripts(String expected, String actual) {
        def expectedLines = normalizeScript(expected)
        def actualLines = normalizeScript(actual)
        
        if (expectedLines.size() != actualLines.size()) {
            return false
        }
        
        // Check if all lines are the same
        for (int i = 0; i < expectedLines.size(); i++) {
            if (expectedLines[i] != actualLines[i]) {
                return false
            }
        }
        
        return true
    }
    
    /**
     * Normalize a script by removing leading/trailing whitespace and empty lines
     */
    List<String> normalizeScript(String script) {
        return script.stripIndent().trim().split('\n').collect { it.trim() }.findAll { it }
    }

    def 'should return base script'() {
        given:
        def expected = '''
            # custom env variables used for azcopy opts
            export AZCOPY_BLOCK_SIZE_MB=4
            export AZCOPY_BLOCK_BLOB_TIER=None
            
            nxf_az_upload() {
                local name=$1
                local target=${2%/} ## remove ending slash
                local base_name="$(basename "$name")"
                local dir_name="$(dirname "$name")"

                local auth_args=""
                if [[ ! -z "$AZCOPY_MSI_CLIENT_ID" ]]; then
                    # When using managed identity, no additional args needed
                    auth_args=""
                else
                    # Use SAS token authentication
                    auth_args="?$AZ_SAS"
                fi

                if [[ -d $name ]]; then
                  if [[ "$base_name" == "$name" ]]; then
                    azcopy cp "$name" "$target$auth_args" --recursive --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
                  else
                    azcopy cp "$name" "$target/$dir_name$auth_args" --recursive --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
                  fi
                else
                  azcopy cp "$name" "$target/$name$auth_args" --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
                fi
            }

            nxf_az_download() {
                local source=$1
                local target=$2
                local basedir=$(dirname $2)
                local ret

                local auth_args=""
                if [[ ! -z "$AZCOPY_MSI_CLIENT_ID" ]]; then
                    # When using managed identity, no additional args needed
                    auth_args=""
                else
                    # Use SAS token authentication
                    auth_args="?$AZ_SAS"
                fi

                mkdir -p "$basedir"

                ret=$(azcopy cp "$source$auth_args" "$target" 2>&1) || {
                    ## if fails check if it was trying to download a directory
                    mkdir -p $target
                    azcopy cp "$source/*$auth_args" "$target" --recursive >/dev/null || {
                        rm -rf $target
                        >&2 echo "Unable to download path: $source"
                        exit 1
                    }
                }
            }
            '''
        
        expect:
        compareScripts(expected, AzBashLib.script())
    }

    def 'should return script with config, with default azcopy opts'() {
        given:
        def opts = new AzCopyOpts()
        def expectedScript = '''
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

            local auth_args=""
            if [[ ! -z "$AZCOPY_MSI_CLIENT_ID" ]]; then
                # When using managed identity, no additional args needed
                auth_args=""
            else
                # Use SAS token authentication
                auth_args="?$AZ_SAS"
            fi

            if [[ -d $name ]]; then
              if [[ "$base_name" == "$name" ]]; then
                azcopy cp "$name" "$target$auth_args" --recursive --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
              else
                azcopy cp "$name" "$target/$dir_name$auth_args" --recursive --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
              fi
            else
              azcopy cp "$name" "$target/$name$auth_args" --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
            fi
        }
        
        nxf_az_download() {
            local source=$1
            local target=$2
            local basedir=$(dirname $2)
            local ret

            local auth_args=""
            if [[ ! -z "$AZCOPY_MSI_CLIENT_ID" ]]; then
                # When using managed identity, no additional args needed
                auth_args=""
            else
                # Use SAS token authentication
                auth_args="?$AZ_SAS"
            fi

            mkdir -p "$basedir"
        
            ret=$(azcopy cp "$source$auth_args" "$target" 2>&1) || {
                ## if fails check if it was trying to download a directory
                mkdir -p $target
                azcopy cp "$source/*$auth_args" "$target" --recursive >/dev/null || {
                    rm -rf $target
                    >&2 echo "Unable to download path: $source"
                    exit 1
                }
            }
        }
        '''
        
        expect:
        compareScripts(expectedScript, AzBashLib.script(opts, 10, 20, Duration.of('30s')))
    }

    def 'should return script with config, with custom azcopy opts'() {
        given:
        def opts = new AzCopyOpts([blobTier: "Hot", blockSize: "10"])
        def expectedScript = '''
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

            local auth_args=""
            if [[ ! -z "$AZCOPY_MSI_CLIENT_ID" ]]; then
                # When using managed identity, no additional args needed
                auth_args=""
            else
                # Use SAS token authentication
                auth_args="?$AZ_SAS"
            fi

            if [[ -d $name ]]; then
              if [[ "$base_name" == "$name" ]]; then
                azcopy cp "$name" "$target$auth_args" --recursive --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
              else
                azcopy cp "$name" "$target/$dir_name$auth_args" --recursive --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
              fi
            else
              azcopy cp "$name" "$target/$name$auth_args" --block-blob-tier $AZCOPY_BLOCK_BLOB_TIER --block-size-mb $AZCOPY_BLOCK_SIZE_MB
            fi
        }
        
        nxf_az_download() {
            local source=$1
            local target=$2
            local basedir=$(dirname $2)
            local ret

            local auth_args=""
            if [[ ! -z "$AZCOPY_MSI_CLIENT_ID" ]]; then
                # When using managed identity, no additional args needed
                auth_args=""
            else
                # Use SAS token authentication
                auth_args="?$AZ_SAS"
            fi

            mkdir -p "$basedir"
        
            ret=$(azcopy cp "$source$auth_args" "$target" 2>&1) || {
                ## if fails check if it was trying to download a directory
                mkdir -p $target
                azcopy cp "$source/*$auth_args" "$target" --recursive >/dev/null || {
                    rm -rf $target
                    >&2 echo "Unable to download path: $source"
                    exit 1
                }
            }
        }
        '''
        
        expect:
        compareScripts(expectedScript, AzBashLib.script(opts, 10, 20, Duration.of('30s')))
    }
}
