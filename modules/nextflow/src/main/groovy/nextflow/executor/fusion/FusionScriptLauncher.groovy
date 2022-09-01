/*
 * Copyright 2020-2022, Seqera Labs
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
 *
 */

package nextflow.executor.fusion

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.BashWrapperBuilder
import nextflow.extension.FilesEx
import nextflow.io.BucketParser
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun
import nextflow.util.Escape
/**
 * Command script launcher implementing the support for Fusion files
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class FusionScriptLauncher extends BashWrapperBuilder {

    private String scheme
    private Path remoteWorkDir
    private Set<String> buckets = new HashSet<>()

    /* ONLY FOR TESTING - DO NOT USE */
    protected FusionScriptLauncher() {
        this.buckets = new HashSet<>()
    }

    FusionScriptLauncher(TaskBean bean, String scheme) {
        super(bean)
        // keep track the google storage work dir
        this.scheme = scheme
        this.remoteWorkDir = bean.workDir

        // map bean work and target dirs to container mount
        // this needed to create the command launcher using container local file paths
        bean.workDir = toContainerMount(bean.workDir)
        bean.targetDir = toContainerMount(bean.targetDir)

        // remap input files to container mounted paths
        for( Map.Entry<String,Path> entry : new HashMap<>(bean.inputFiles).entrySet() ) {
            bean.inputFiles.put( entry.key, toContainerMount(entry.value) )
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
        return "NXF_CHDIR=${Escape.path(bean.workDir)}\n"
    }

    Path toContainerMount(Path path) {
        if( path == null )
            return null

        final p = BucketParser.from( FilesEx.toUriString(path) )

        if( p.scheme != scheme )
            throw new IllegalArgumentException("Unexpected path for Fusion script launcher: ${path.toUriString()}")

        final result = "/fusion/$p.scheme/${p.bucket}${p.path}"
        buckets.add(p.bucket)
        return Path.of(result)
    }

    Set<String> fusionBuckets() {
        return buckets
    }

    String fusionWork() {
        return toContainerMount(remoteWorkDir).toString()
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
}
