package nextflow.executor

import java.nio.file.Paths

import nextflow.Global
import nextflow.Session
import nextflow.cloud.azure.batch.AzBatchExecutor
import nextflow.cloud.azure.batch.AzFileCopyStrategy
import nextflow.cloud.azure.batch.AzHelper
import nextflow.cloud.azure.config.AzConfig
import nextflow.cloud.azure.file.AzPathFactory
import nextflow.processor.TaskBean
import spock.lang.Requires
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BashWrapperBuilderWithAzTest extends Specification {

    @Requires({System.getenv('AZURE_STORAGE_ACCOUNT_NAME') && System.getenv('AZURE_STORAGE_ACCOUNT_KEY')})
    def 'should include az helpers' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [:] }
        and:
        def folder = Paths.get('/work/dir')
        def target = AzPathFactory.parse('az://some/bucket')

        def bean = new TaskBean([
                name: 'Hello 1',
                workDir: folder,
                targetDir: target,
                scratch: true,
                outputFiles: ['test.bam','test.bai'],
                script: 'echo Hello world!',
        ])

        def copy = new SimpleFileCopyStrategy(bean)

        /*
         * simple bash run
         */
        when:
        def binding = new BashWrapperBuilder(bean,copy).makeBinding()
        then:
        binding.unstage_outputs == """\
                    IFS=\$'\\n'
                    for name in \$(eval "ls -1d test.bam test.bai" | sort | uniq); do
                        nxf_az_upload \$name '${AzHelper.toHttpUrl(target)}' || true
                    done
                    unset IFS
                    """.stripIndent().rightTrim()

        binding.helpers_script == '''\
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
            
            '''.stripIndent(true)
    }

    @Requires({System.getenv('AZURE_STORAGE_ACCOUNT_NAME') && System.getenv('AZURE_STORAGE_ACCOUNT_KEY')})
    def 'should include az helpers and bash lib' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [:] }
        and:
        def folder = AzPathFactory.parse('az://work/dir')
        def target = AzPathFactory.parse('az://some/bucket')

        def bean = new TaskBean([
                name: 'Hello 1',
                workDir: folder,
                targetDir: target,
                scratch: true,
                outputFiles: ['test.bam','test.bai'],
                script: 'echo Hello world!',
        ])

        def exec = Mock(AzBatchExecutor) {
            getConfig() >> new AzConfig([:])
        }
        and:
        def copy = new AzFileCopyStrategy(bean, exec)

        /*
         * simple bash run
         */
        when:
        def binding = new BashWrapperBuilder(bean,copy).makeBinding()
        then:
        binding.unstage_outputs == """\
                    uploads=()
                    IFS=\$'\\n'
                    for name in \$(eval "ls -1d test.bam test.bai" | sort | uniq); do
                        uploads+=("nxf_az_upload '\$name' '${AzHelper.toHttpUrl(target)}'")
                    done
                    unset IFS
                    nxf_parallel "\${uploads[@]}"
                    """.stripIndent()

        binding.helpers_script == '''\
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
            
            '''.stripIndent(true)
    }

}
