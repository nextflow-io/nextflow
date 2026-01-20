/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.cli

import groovy.transform.CompileStatic
import org.pf4j.ExtensionPoint

/**
 * Extension point interface for the `launch` command.
 *
 * @see io.seqera.tower.plugin.launch.LaunchCommandImpl
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
interface LaunchCommand extends ExtensionPoint {
    void launch(LaunchOptions options)
}


/**
 * Data class to hold launch options
 */
@CompileStatic
class LaunchOptions {
    String pipeline
    String workspace
    String computeEnv
    String runName
    String workDir
    String revision
    String profile
    List<String> configFiles
    String paramsFile
    String entryName
    String resume
    boolean latest
    boolean stubRun
    String mainScript
    Map<String, String> params
}
