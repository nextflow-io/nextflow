/*
 * Copyright 2023, Seqera Labs
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

package nextflow.cloud.google.batch


import com.google.cloud.batch.v1.Volume
import groovy.transform.CompileStatic
import nextflow.fusion.FusionAwareTask
import nextflow.fusion.FusionScriptLauncher
/**
 * A mere adapter for {@link GoogleBatchLauncherSpec} for fusion
 * launch
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class GoogleBatchFusionAdapter implements GoogleBatchLauncherSpec {

    private FusionAwareTask task

    private FusionScriptLauncher launcher

    GoogleBatchFusionAdapter(FusionAwareTask task, FusionScriptLauncher launcher) {
        this.task = task
        this.launcher = launcher
    }

    @Override
    List<String> getContainerMounts() {
        return List.of()
    }

    @Override
    List<Volume> getVolumes() {
        return [
                Volume.newBuilder()
                        .setDeviceName('fusion')
                        .setMountPath('/tmp')
                        .build()
        ]
    }

    @Override
    String runCommand() {
        throw new UnsupportedOperationException()
    }

    @Override
    List<String> launchCommand() {
        return task.fusionSubmitCli()
    }

    @Override
    Map<String, String> getEnvironment() {
        return launcher.fusionEnv()
    }
}
