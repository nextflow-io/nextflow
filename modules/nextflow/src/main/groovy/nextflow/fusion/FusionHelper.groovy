/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.fusion

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.Session
import nextflow.SysEnv
import nextflow.container.ContainerBuilder
import nextflow.container.ContainerConfig
import nextflow.extension.FilesEx
import nextflow.io.BucketParser

/**
 * Helper method to handle fusion common logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class FusionHelper {

    @Memoized
    static boolean isFusionEnabled(Session session) {
        def result = session.config.navigate('fusion.enabled')
        if( result == null )
            result = SysEnv.get('FUSION_ENABLED')
        return result!=null ? result.toString()=='true' : false
    }

    static String runWithContainer(FusionScriptLauncher launcher, ContainerConfig containerConfig, String containerName, String containerOpts, List<String> runCmd) {
        if( !containerName )
            throw new IllegalArgumentException("Missing task container -- Fusion requires the task to be executed by a container process")
        final engine = containerConfig.getEngine()
        final containerBuilder = ContainerBuilder.create(engine, containerName)
                .addMountWorkDir(false)
                .addRunOptions(containerConfig.runOptions as String)
                .addRunOptions(containerOpts)
                .addRunOptions(containerConfig.fusionOptions())
                .params(containerConfig)

        // add fusion env vars
        for(Map.Entry<String,String> it : launcher.fusionEnv()) {
            containerBuilder.addEnv("$it.key=$it.value")
        }

        // add env variables
        for( String env : containerConfig.getEnvWhitelist())
            containerBuilder.addEnv(env)

        // patch the cmd wrapping the last item in quotes
        final patchCmd = new ArrayList(runCmd)
        patchCmd[-1] = "'${patchCmd[-1]}'"

        // assemble the final command
        final containerCmd = containerBuilder
                .build()
                .getRunCommand(patchCmd.join(' '))
                .replaceAll('-w "\\$NXF_TASK_WORKDIR" ','') // <-- hack to remove the PWD work dir

        return containerCmd
    }

    static Path toContainerMount(Path path, String scheme) {
        if( path == null )
            return null

        final p = BucketParser.from( FilesEx.toUriString(path) )

        if( p.scheme != scheme )
            throw new IllegalArgumentException("Unexpected path for Fusion script launcher: ${path.toUriString()}")

        final result = "/fusion/$p.scheme/${p.bucket}${p.path}"
        return Path.of(result)
    }

    static Path toContainerMount(Path path) {
        toContainerMount(path, path.scheme)
    }

}
