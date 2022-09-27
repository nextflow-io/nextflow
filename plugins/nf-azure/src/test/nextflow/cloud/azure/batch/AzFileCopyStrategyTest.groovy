package nextflow.cloud.azure.batch

import java.nio.file.FileSystem
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.spi.FileSystemProvider

import com.azure.storage.blob.BlobClient
import nextflow.Session
import nextflow.cloud.azure.config.AzConfig
import nextflow.cloud.azure.nio.AzPath
import nextflow.processor.TaskBean
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzFileCopyStrategyTest extends Specification {

    protected Path mockAzPath(String path, boolean isDir=false) {
        assert path.startsWith('az://')

        def tokens = path.tokenize('/')
        def bucket = tokens[1]
        def file = '/' + tokens[2..-1].join('/')

        def attr = Mock(BasicFileAttributes)
        attr.isDirectory() >> isDir
        attr.isRegularFile() >> !isDir
        attr.isSymbolicLink() >> false

        def provider = Mock(FileSystemProvider)
        provider.getScheme() >> 'az'
        provider.readAttributes(_, _, _) >> attr

        def fs = Mock(FileSystem)
        fs.provider() >> provider
        fs.toString() >> ('az://' + bucket)
        def uri = GroovyMock(URI)
        uri.toString() >> path

        def BLOB_CLIENT = Mock(BlobClient) {
            getBlobUrl() >> { path.replace('az://', 'http://account.blob.core.windows.net/') }
        }

        def result = GroovyMock(AzPath)
        result.toUriString() >> path
        result.toString() >> file
        result.getFileSystem() >> fs
        result.toUri() >> uri
        result.resolve(_) >> { mockAzPath("$path/${it[0]}") }
        result.toAbsolutePath() >> result
        result.asBoolean() >> true
        result.getParent() >> { def p=path.lastIndexOf('/'); p!=-1 ? mockAzPath("${path.substring(0,p)}", true) : null }
        result.getFileName() >> { Paths.get(tokens[-1]) }
        result.getName() >> tokens[1]
        result.blobClient() >> BLOB_CLIENT
        return result
    }


    def setup() {
        new Session()
    }

    def 'test bash wrapper with input'() {
        given:
        def workDir = mockAzPath( 'az://my-data/work/dir' )
        def token = '12345'
        def config = new AzConfig([storage:[sasToken: token]])
        def executor = Mock(AzBatchExecutor) { getConfig() >> config }

        when:
        def binding = new AzBatchScriptLauncher([
                name: 'Hello 1',
                workDir: workDir,
                script: 'echo Hello world!',
                environment: [FOO: 1, BAR:'any'],
                input: 'Ciao ciao' ] as TaskBean, executor) .makeBinding()

        then:
        binding.stage_inputs == '''\
                # stage input files
                downloads=(true)
                
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()

        binding.unstage_controls == '''\
                nxf_az_upload .command.out 'http://account.blob.core.windows.net/my-data/work/dir' || true
                nxf_az_upload .command.err 'http://account.blob.core.windows.net/my-data/work/dir' || true
                '''.stripIndent()

        binding.launch_cmd == '/bin/bash -ue .command.sh < .command.in'
        binding.unstage_outputs == null

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

    def 'should include remote bind dir' () {
        given:
        def remoteBin = mockAzPath( 'az://my-data/work/remote/bin' )
        def workDir = mockAzPath( 'az://my-data/work/dir' )
        def token = '12345'
        def config = new AzConfig([storage:[sasToken: token]])
        def executor = Mock(AzBatchExecutor) {
            getConfig() >> config
            getRemoteBinDir() >> remoteBin
        }

        when:
        def binding = new AzBatchScriptLauncher([
                name: 'Hello 1',
                workDir: workDir,
                script: 'echo Hello world!',
                environment: [FOO: 1, BAR:'any'],
                input: 'Ciao ciao' ] as TaskBean, executor) .makeBinding()

        then:
        binding.stage_inputs == '''\
                # stage input files
                nxf_az_download 'http://account.blob.core.windows.net/my-data/work/remote/bin' $PWD/.nextflow-bin
                chmod +x $PWD/.nextflow-bin/* || true
                downloads=(true)
                
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()

        binding.task_env == '''\
                    export FOO="1"
                    export BAR="any"
                    export PATH="$PWD/.nextflow-bin:$AZ_BATCH_NODE_SHARED_DIR/bin/:$PATH"
                    export AZCOPY_LOG_LOCATION="$PWD/.azcopy_log"
                    export AZ_SAS="12345"
                    '''.stripIndent()

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



    def 'test bash wrapper with outputs and stats'() {
        /*
         * simple bash run
         */
        given:
        def workDir = mockAzPath( 'az://my-data/work/dir' )
        def input1 = mockAzPath('az://my-data/work/dir/file1.txt')
        def input2 = mockAzPath('az://my-data/work/dir/file2.txt')
        def token = '12345'
        def config = new AzConfig([storage:[sasToken: token]])
        def executor = Mock(AzBatchExecutor) { getConfig() >> config }

        when:
        def binding = new AzBatchScriptLauncher([
                name: 'Hello 1',
                workDir: workDir,
                targetDir: workDir,
                statsEnabled: true,
                inputFiles: ['file1.txt': input1, 'file2.txt': input2],
                outputFiles: ['foo.txt', 'bar.fastq'],
                script: 'echo Hello world!',
                input: 'Ciao ciao' ] as TaskBean, executor) .makeBinding()

        then:

        binding.unstage_controls == '''\
                nxf_az_upload .command.out 'http://account.blob.core.windows.net/my-data/work/dir' || true
                nxf_az_upload .command.err 'http://account.blob.core.windows.net/my-data/work/dir' || true
                nxf_az_upload .command.trace 'http://account.blob.core.windows.net/my-data/work/dir' || true
                '''.stripIndent()

        binding.stage_inputs == '''\
                # stage input files
                downloads=(true)
                rm -f file1.txt
                rm -f file2.txt
                downloads+=("nxf_az_download 'http://account.blob.core.windows.net/my-data/work/dir/file1.txt' file1.txt")
                downloads+=("nxf_az_download 'http://account.blob.core.windows.net/my-data/work/dir/file2.txt' file2.txt")
                nxf_parallel "${downloads[@]}"
                '''.stripIndent()

        binding.unstage_outputs == '''
                    uploads=()
                    IFS=$'\\n'
                    for name in $(eval "ls -1d foo.txt bar.fastq" | sort | uniq); do
                        uploads+=("nxf_az_upload '$name' 'http://account.blob.core.windows.net/my-data/work/dir'")
                    done
                    unset IFS
                    nxf_parallel "${uploads[@]}"
                    '''.stripIndent().leftTrim()

        binding.launch_cmd == '/bin/bash .command.run nxf_trace'

        binding.task_env == '''\
                    export PATH="$PWD/.nextflow-bin:$AZ_BATCH_NODE_SHARED_DIR/bin/:$PATH"
                    export AZCOPY_LOG_LOCATION="$PWD/.azcopy_log"
                    export AZ_SAS="12345"
                    '''.stripIndent()

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
