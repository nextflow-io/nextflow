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

package nextflow.cloud.google.batch

import java.nio.file.Path

import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.google.batch.client.BatchConfig
import nextflow.executor.SimpleFileCopyStrategy
import nextflow.extension.FilesEx
import nextflow.processor.TaskBean
import nextflow.util.Escape

/**
 * Optional Google Batch staging when {@link BatchConfig#usesGoogleBatchStaging()} is true, honouring
 * {@link BatchConfig#stageInCopyTransport} / {@link BatchConfig#stageOutCopyTransport} for {@code copy} modes.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class GoogleBatchFileCopyStrategy extends SimpleFileCopyStrategy {

    private final BatchConfig batchConfig

    GoogleBatchFileCopyStrategy(TaskBean bean, BatchConfig batchConfig) {
        super(bean)
        this.batchConfig = batchConfig
    }

    /**
     * CLI stage-in when a CLI transport is selected and {@code stageInMode=copy}.
     */
    private boolean useCliStageInCopy() {
        return BatchConfig.isCliCopyTransport(batchConfig.stageInCopyTransport) && 'copy' == stageinMode
    }

    private String effectiveStageOutMode() {
        stageoutMode ?: ( workDir==targetDir ? 'copy' : 'move' )
    }

    private boolean shouldUseCliStageOutCopy() {
        if( effectiveStageOutMode() != 'copy' )
            return false
        return BatchConfig.isCliCopyTransport(batchConfig.stageOutCopyTransport)
    }

    private boolean needsGoogleBatchBashLib() {
        useCliStageInCopy() || shouldUseCliStageOutCopy()
    }

    @Override
    String getBeforeStartScript() {
        final base = super.getBeforeStartScript()
        if( !needsGoogleBatchBashLib() )
            return base
        final gs = GoogleBatchBashLib.script(batchConfig)
        return gs + (base ? '\n' + base : '')
    }

    @Override
    String getStageInputFilesScript(Map<String,Path> inputFiles) {
        if( !useCliStageInCopy() )
            return super.getStageInputFilesScript(inputFiles)
        def result = 'downloads=(true)\n'
        result += super.getStageInputFilesScript(inputFiles) + '\n'
        result += 'nxf_parallel "${downloads[@]}"\n'
        return result
    }

    @Override
    String stageInputFile(Path path, String targetName) {
        if( !useCliStageInCopy() )
            return super.stageInputFile(path, targetName)
        final gsUri = gsUriForCliStageIn(path)
        if( gsUri != null ) {
            final cmd = batchConfig.maxTransferAttempts > 1
                    ? "downloads+=(\"nxf_cp_retry nxf_gs_download '${gsUri}' ${Escape.path(targetName)}\")"
                    : "downloads+=(\"nxf_gs_download '${gsUri}' ${Escape.path(targetName)}\")"
            return cmd
        }
        return super.stageInputFile(path, targetName)
    }

    /**
     * Resolve a {@code gs://} URI for CLI stage-in: from {@link CloudStoragePath}, or from a container fuse path under {@link GoogleBatchScriptLauncher#MOUNT_ROOT}.
     */
    private static String gsUriForCliStageIn(Path path) {
        if( path instanceof CloudStoragePath )
            return FilesEx.toUriString((CloudStoragePath)path)
        final s = path.toString()
        final prefix = GoogleBatchScriptLauncher.MOUNT_ROOT + '/'
        if( s.startsWith(prefix) )
            return toGsUriFromContainerMount(path)
        return null
    }

    @Override
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir) {
        final mode = effectiveStageOutMode()
        if( mode == 'move' )
            return super.getUnstageOutputFilesScript(outputFiles, targetDir)
        if( mode == 'copy' && shouldUseCliStageOutCopy() )
            return getUnstageOutputFilesScriptGcloud(outputFiles, targetDir)
        return super.getUnstageOutputFilesScript(outputFiles, targetDir)
    }

    private String getUnstageOutputFilesScriptGcloud(List<String> outputFiles, Path targetDir) {
        final patterns = normalizeGlobStarPaths(outputFiles)
        log.trace "[GOOGLE BATCH] Unstaging file path (CLI transport): $patterns"

        if( !patterns )
            return null

        final gsTarget = toGsUriFromContainerMount(targetDir)
        final escape = new ArrayList(outputFiles.size())
        for( String it : patterns )
            escape.add( Escape.path(it) )

        return """\
            uploads=()
            IFS=\$'\\n'
            for name in \$(eval "ls -1d ${escape.join(' ')}" | sort | uniq); do
                uploads+=("nxf_gs_upload '\$name' '${gsTarget}' 0")
            done
            unset IFS
            nxf_parallel "\${uploads[@]}"
            """.stripIndent(true)
    }

    @Override
    String fileStr(Path path) {
        !useCliStageInCopy() ? super.fileStr(path) : Escape.path(path.getFileName())
    }

    @Override
    String pipeInputFile(Path file) {
        !useCliStageInCopy() ? super.pipeInputFile(file) : " < ${Escape.path(file.getFileName())}"
    }

    static String toGsUriFromContainerMount(Path containerPath) {
        final s = containerPath.toString()
        final prefix = GoogleBatchScriptLauncher.MOUNT_ROOT + '/'
        if( !s.startsWith(prefix) )
            throw new IllegalArgumentException("Expected path under ${GoogleBatchScriptLauncher.MOUNT_ROOT}, got: $s")
        final rest = s.substring(prefix.length())
        final slash = rest.indexOf('/')
        if( slash < 0 )
            return "gs://$rest/"
        final bucket = rest.substring(0, slash)
        final obj = rest.substring(slash)
        return "gs://$bucket$obj"
    }
}
