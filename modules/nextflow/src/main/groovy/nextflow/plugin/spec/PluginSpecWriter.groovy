/*
 * Copyright 2024-2025, Seqera Labs
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
package nextflow.plugin.spec

import groovy.json.JsonOutput
import groovy.transform.CompileStatic

/**
 * This entrypoint is used by the Nextflow Gradle plugin
 * to generate plugin specs when packaging a plugin.
 *
 * It must be defined in the Nextflow runtime so that the
 * Gradle plugin can use it without depending on Nextflow at
 * compile-time (i.e. the plugin should use the version of
 * Nextflow specified by the user).
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class PluginSpecWriter {

    static void main(String[] args) {
        if (args.length < 2) {
            System.err.println("Usage: PluginSpecWriter <output-path> <class1> [class2] ...")
            System.exit(1)
        }

        final outputPath = args[0]
        final extensionPoints = args[1..-1]

        // build plugin spec
        final spec = new PluginSpec(extensionPoints).build()

        // write plugin spec to JSON file
        final file = new File(outputPath)
        file.parentFile.mkdirs()
        file.text = JsonOutput.toJson(spec)

        println "Saved plugin spec to $file"
    }
}
