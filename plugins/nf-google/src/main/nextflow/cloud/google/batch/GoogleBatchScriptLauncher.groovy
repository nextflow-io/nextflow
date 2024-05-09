/*
 * Copyright 2023, Seqera Labs.
 * Copyright 2022, Google Inc.
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
import java.nio.file.Paths

import com.google.cloud.batch.v1.GCS
import com.google.cloud.batch.v1.Volume
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.google.batch.client.BatchConfig
import nextflow.executor.BashWrapperBuilder
import nextflow.extension.FilesEx
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun
import nextflow.util.Escape
import nextflow.util.PathTrie

/**
 * Implement Nextflow task launcher script
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class GoogleBatchScriptLauncher extends BashWrapperBuilder implements GoogleBatchLauncherSpec {

    private static final String MOUNT_ROOT = '/mnt/disks'

    private BatchConfig config
    private CloudStoragePath remoteWorkDir
    private Path remoteBinDir
    private Set<String> buckets = new HashSet<>()
    private PathTrie pathTrie = new PathTrie()

    /* ONLY FOR TESTING - DO NOT USE */
    protected GoogleBatchScriptLauncher() {}

    GoogleBatchScriptLauncher(TaskBean bean, Path remoteBinDir) {
        super(bean)
        // keep track the google storage work dir
        this.remoteWorkDir = (CloudStoragePath) bean.workDir
        this.remoteBinDir = toContainerMount(remoteBinDir)

        // map bean work and target dirs to container mount
        // this is needed to create the command launcher using container local file paths
        bean.workDir = toContainerMount(bean.workDir)
        bean.targetDir = toContainerMount(bean.targetDir)

        // add all children work dir 
        if( bean.arrayWorkDirs ) {
            for( Path it : bean.arrayWorkDirs )
                toContainerMount(it)
        }

        // remap input files to container mounted paths
        for( Map.Entry<String,Path> entry : new HashMap<>(bean.inputFiles).entrySet() ) {
            bean.inputFiles.put( entry.key, toContainerMount(entry.value, true) )
        }

        // include task script as an input to force its staging in the container work directory
        bean.inputFiles[TaskRun.CMD_SCRIPT] = bean.workDir.resolve(TaskRun.CMD_SCRIPT)
        // add the wrapper file when stats are enabled
        // NOTE: this must match the logic that uses the run script in BashWrapperBuilder
        if( isTraceRequired() ) {
            bean.inputFiles[TaskRun.CMD_RUN] = bean.workDir.resolve(TaskRun.CMD_RUN)
        }
        // include task stdin file
        if( bean.input != null ) {
            bean.inputFiles[TaskRun.CMD_INFILE] = bean.workDir.resolve(TaskRun.CMD_INFILE)
        }

        // make it change to the task work dir
        bean.headerScript = headerScript(bean)
        // enable use of local scratch dir
        if( scratch==null )
            scratch = true
    }

    protected String headerScript(TaskBean bean) {
        def result = "NXF_CHDIR=${Escape.path(bean.workDir)}\n"
        if( remoteBinDir ) {
            result += "cp -r $remoteBinDir \$HOME/.nextflow-bin\n"
            result += 'chmod +x $HOME/.nextflow-bin/*\n'
            result += 'export PATH=$HOME/.nextflow-bin:$PATH\n'
        }
        return result
    }

    protected Path toContainerMount(Path path, boolean parent=false) {
        if( path instanceof CloudStoragePath ) {
            buckets.add(path.bucket())
            pathTrie.add( (parent ? "/${path.bucket()}${path.parent}" : "/${path.bucket()}${path}").toString() )
            final containerMount = containerMountPath(path)
            log.trace "Path ${FilesEx.toUriString(path)} to container mount: $containerMount"
            return Paths.get(containerMount)
        }
        else if( path==null )
            return null
        throw new IllegalArgumentException("Unexpected path for Google Batch task handler: ${path.toUriString()}")
    }

    @Override
    String runCommand() {
        launchCommand(workDirMount)
    }

    @Override
    List<String> getContainerMounts() {
        final result = new ArrayList(10)
        for( String it : pathTrie.longest() ) {
            result.add( "${MOUNT_ROOT}${it}:${MOUNT_ROOT}${it}:rw".toString() )
        }
        return result
    }

    @Override
    List<Volume> getVolumes() {
        final result = new ArrayList(10)
        for( String it : buckets ) {
            final mountOptions = ['-o rw', '-implicit-dirs']
            if( config && config.googleOpts.enableRequesterPaysBuckets )
                mountOptions << "--billing-project ${config.googleOpts.projectId}".toString()

            result.add(
                Volume.newBuilder()
                    .setGcs(
                        GCS.newBuilder()
                            .setRemotePath(it)
                    )
                    .setMountPath( "${MOUNT_ROOT}/${it}".toString() )
                    .addAllMountOptions( mountOptions )
                    .build()
            )
        }
        return result
    }

    String getWorkDirMount() {
        return workDir.toString()
    }

    @Override
    protected Path targetWrapperFile() {
        return remoteWorkDir.resolve(TaskRun.CMD_RUN)
    }

    @Override
    protected Path targetScriptFile() {
        return remoteWorkDir.resolve(TaskRun.CMD_SCRIPT)
    }

    @Override
    protected Path targetInputFile() {
        return remoteWorkDir.resolve(TaskRun.CMD_INFILE)
    }

    GoogleBatchScriptLauncher withConfig(BatchConfig config) {
        this.config = config
        return this
    }

    static String launchCommand( String workDir ) {
        "trap \"{ cp ${TaskRun.CMD_LOG} ${workDir}/${TaskRun.CMD_LOG}; }\" ERR; /bin/bash ${workDir}/${TaskRun.CMD_RUN} 2>&1 | tee ${TaskRun.CMD_LOG}"
    }

    static String containerMountPath(CloudStoragePath path) {
        return "$MOUNT_ROOT/${path.bucket()}${path}"
    }
}
