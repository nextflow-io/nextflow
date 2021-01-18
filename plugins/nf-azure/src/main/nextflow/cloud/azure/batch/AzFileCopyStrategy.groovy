package nextflow.cloud.azure.batch

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.CloudTransferOptions
import nextflow.executor.BashFunLib
import nextflow.executor.SimpleFileCopyStrategy
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.Escape
/**
 * Implements file copy strategy for Azure Batch
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AzFileCopyStrategy extends SimpleFileCopyStrategy {

    private int maxTransferAttempts
    private int maxParallelTransfers
    private Duration delayBetweenAttempts
    private String rcloneCli = './rclone'
    private String sas

    protected AzFileCopyStrategy() {}

    AzFileCopyStrategy(TaskBean bean, String sas) {
        super(bean)
        this.sas = sas
        //
        this.maxParallelTransfers = 1 // /*config.maxParallelTransfers ?:*/ CloudTransferOptions.MAX_TRANSFER
        this.maxTransferAttempts = /*config.maxTransferAttempts ?:*/ CloudTransferOptions.MAX_TRANSFER_ATTEMPTS
        this.delayBetweenAttempts = /*config.delayBetweenAttempts ?:*/ CloudTransferOptions.DEFAULT_DELAY_BETWEEN_ATTEMPTS
    }

    protected String toAzHttpUrl(Path path) {
        AzHelper.toHttpUrl(path, sas)
    }

    @Override
    String getBeforeStartScript() {

        BashFunLib.body(maxParallelTransfers, maxTransferAttempts, delayBetweenAttempts) +
                """
        SAS='$sas'

        nxf_az_upload() {
            local name=\$1
            local target=\${2%/} ## remove ending slash
        
            if [[ -d \$name ]]; then
              $rcloneCli copy "\$name" ":azureblob:/\$target" --azureblob-sas-url \"\$SAS\"
            else 
              $rcloneCli copy "\$name" ":azureblob:/\$target/\$name --azureblob-sas-url \"\$SAS\""
            fi  
        }
        
        nxf_az_download() {
            local source=\$1
            local target=\$2
            local basedir=\$(dirname \$2)
            local ret
            mkdir -p "\$basedir"
        
            ret=\$($rcloneCli copy :azureblob\$source "\$target" --azureblob-sas-url \"\$SAS\" 2>&1) || {
                ## if fails check if it was trying to download a directory
                mkdir \$target
                $rcloneCli copy ":azureblob:\$source/*" "\$target" --azureblob-sas-url \"\$SAS\" >/dev/null || {
                    rm -rf \$target
                    >&2 echo "Unable to download path: \$source"
                    exit 1
                }
            }
        }
        """.stripIndent()
    }

    @Override
    String getStageInputFilesScript(Map<String, Path> inputFiles) {
        def result = 'downloads=()\n'
        result += super.getStageInputFilesScript(inputFiles) + '\n'
        result += 'nxf_parallel "${downloads[@]}"\n'
        return result
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String stageInputFile( Path path, String targetName ) {
        // third param should not be escaped, because it's used in the grep match rule
        def stage_cmd = maxTransferAttempts > 1
                ? "downloads+=(\"nxf_cp_retry nxf_az_download '${toAzHttpUrl(path)}' ${Escape.path(targetName)}\")"
                : "downloads+=(\"nxf_az_download '${toAzHttpUrl(path)}' ${Escape.path(targetName)}\")"
        return stage_cmd
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir) {
        final patterns = normalizeGlobStarPaths(outputFiles)
        // create a bash script that will copy the out file to the working directory
        log.trace "[AZURE BATCH] Unstaging file path: $patterns"

        if( !patterns )
            return null

        final escape = new ArrayList(outputFiles.size())
        for( String it : patterns )
            escape.add( Escape.path(it) )

        return """\
            uploads=()
            IFS=\$'\\n'
            for name in \$(eval "ls -1d ${escape.join(' ')}" | sort | uniq); do
                uploads+=("nxf_az_upload '\$name' '${toAzHttpUrl(targetDir)}'")
            done
            unset IFS
            nxf_parallel "\${uploads[@]}"
            """.stripIndent(true)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String touchFile( Path file ) {
        "echo start > .command.begin && nxf_az_upload .command.begin '${toAzHttpUrl(file.parent)}'"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String fileStr( Path path ) {
        Escape.path(path.getFileName())
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String copyFile( String name, Path target ) {
        "nxf_az_upload ${Escape.path(name)} '${toAzHttpUrl(target.parent)}'"
    }

    /**
     * {@inheritDoc}
     */
    String exitFile( Path path ) {
        " > ${TaskRun.CMD_EXIT} && nxf_az_upload ${TaskRun.CMD_EXIT} '${toAzHttpUrl(path.parent)}' || true"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String pipeInputFile( Path path ) {
        " < ${Escape.path(path.getFileName())}"
    }

}
