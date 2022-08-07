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
import groovy.util.logging.Slf4j
import org.pf4j.BasePluginLoader
import org.pf4j.PluginManager

/**
 * Plugin loader specialised for unit testing. The plugin must be defined
 * in the sub-project 'testFixtures' source tree
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TestPluginLoader extends BasePluginLoader {

    TestPluginLoader(PluginManager pluginManager) {
        super(pluginManager, new TestPluginClasspath())
    }

}
