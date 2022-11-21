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

import static nextflow.executor.fusion.FusionHelper.*

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.executor.BashWrapperBuilder
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
    private Map<String,String> env

    /* ONLY FOR TESTING - DO NOT USE */
    protected FusionScriptLauncher() {
        this.buckets = new HashSet<>()
    }

    static FusionScriptLauncher create(TaskBean bean, String scheme) {

        final buckets = new HashSet(10)
        final remoteWorkDir = bean.workDir

        // map bean work and target dirs to container mount
        // this needed to create the command launcher using container local file paths
        bean.workDir = toContainerMount(bean.workDir, scheme, buckets)
        bean.targetDir = toContainerMount(bean.targetDir, scheme, buckets)

        // remap input files to container mounted paths
        for( Map.Entry<String,Path> entry : new HashMap<>(bean.inputFiles).entrySet() ) {
            bean.inputFiles.put( entry.key, toContainerMount(entry.value, scheme, buckets) )
        }

        // make it change to the task work dir
        bean.headerScript = headerScript(bean)
        // enable use of local scratch dir
        if( bean.scratch==null )
            bean.scratch = true

        return new FusionScriptLauncher(bean, scheme, remoteWorkDir, buckets)
    }

    FusionScriptLauncher(TaskBean bean, String scheme, Path remoteWorkDir, Set<String> buckets) {
        super(bean)
        // keep track the google storage work dir
        this.scheme = scheme
        this.remoteWorkDir = remoteWorkDir
        this.buckets = buckets
    }

    static protected String headerScript(TaskBean bean) {
        return "NXF_CHDIR=${Escape.path(bean.workDir)}\n"
    }

    Path toContainerMount(Path path) {
        toContainerMount(path,scheme,buckets)
    }

    Set<String> fusionBuckets() {
        return buckets
    }

    Map<String,String> fusionEnv() {
        if( env==null ) {
            final buckets = fusionBuckets().collect(it->"$scheme://$it").join(',')
            final work = toContainerMount(remoteWorkDir).toString()
            final result = new LinkedHashMap(10)
            result.NXF_FUSION_WORK = work
            result.NXF_FUSION_BUCKETS = buckets
            final endpoint = Global.getAwsS3Endpoint()
            final creds = exportAwsAccessKeys() ? Global.getAwsCredentials() : Collections.<String>emptyList()
            if( creds ) {
                result.AWS_ACCESS_KEY_ID = creds[0]
                result.AWS_SECRET_ACCESS_KEY = creds[1]
            }
            if( endpoint )
                result.AWS_S3_ENDPOINT = endpoint
            env = result
        }
        return env
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

    boolean exportAwsAccessKeys() {
        exportAwsAccessKeys0()
    }

    @Memoized
    protected boolean exportAwsAccessKeys0() {
        return Global.config?.navigate('fusion.exportAwsAccessKeys', false)
    }
}
