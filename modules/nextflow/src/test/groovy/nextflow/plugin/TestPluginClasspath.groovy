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

package nextflow.plugin

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.pf4j.PluginClasspath

/**
 * A plugin classpath that looks for the plugin in the
 * 'testPlugin' source set classes and resources
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TestPluginClasspath extends PluginClasspath {

    private boolean logged

    TestPluginClasspath() {
        // the output dirs of the 'testPlugin' source set, provided by the Gradle test task
        addClassesDirectories(System.getProperty('nextflow.testPlugin.classesDirs').tokenize(File.pathSeparator))
        addClassesDirectories(System.getProperty('nextflow.testPlugin.resourcesDir'))
    }

    @Override
    Set<String> getClassesDirectories() {
        if( !logged ) {
            log.debug "Groovy TEST plugin classpath: classes-dirs=${super.getClassesDirectories()}; jars-dirs=${super.getJarsDirectories()}"
            logged = true
        }
        return super.getClassesDirectories()
    }

    @Override
    Set<String> getJarsDirectories() {
        if( !logged ) {
            log.debug "Groovy TEST plugin classpath: classes-dirs=${super.getClassesDirectories()}; jars-dirs=${super.getJarsDirectories()}"
            logged = true
        }
        return super.getJarsDirectories()
    }
}
