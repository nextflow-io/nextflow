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

package nextflow.plugin

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j;
import org.pf4j.PluginClasspath;

/**
 * Customise classpath loader for Groovy based
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class DevPluginClasspath extends PluginClasspath {

    private boolean logged

    DevPluginClasspath() {
        // the path where classes are resources should be found in the dev environment
        // for each plugin project directory
        addClassesDirectories("build/classes/groovy/main", "build/resources/main", 'build/classes/main')

        // the path where the plugin dependencies jar files are expected to be found
        // note: this path is not created automatically by Gradle, it should be created by a custom task
        // see `targetLibs` task in the base plugins `build.gradle`
        addJarsDirectories('build/target/libs')

    }

    @Override
    Set<String> getClassesDirectories() {
        if( !logged ) {
            log.debug "Groovy DEV plugin classpath: classes-dirs=${super.getClassesDirectories()}; jars-dirs=${super.getJarsDirectories()}"
            logged = true
        }
        return super.getClassesDirectories()
    }

    @Override
    Set<String> getJarsDirectories() {
        if( !logged ) {
            log.debug "Groovy DEV plugin classpath: classes-dirs=${super.getClassesDirectories()}; jars-dirs=${super.getJarsDirectories()}"
            logged = true
        }
        return super.getJarsDirectories()
    }
}
