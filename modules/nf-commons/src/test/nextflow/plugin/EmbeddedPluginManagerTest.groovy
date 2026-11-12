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

import java.nio.file.Files

import spock.lang.Specification

import static nextflow.plugin.LocalPluginManagerTest.createPlugin

class EmbeddedPluginManagerTest extends Specification {

    def 'should load each plugin in namespace folders once' () {
        given:
        def store = Files.createTempDirectory('test')
        createPlugin(store, 'nf-bar', '1.0.0')
        createPlugin(store, 'org1/nf-foo', '1.0.0')
        createPlugin(store, 'org1/nf-baz', '1.0.0')
        createPlugin(store, 'org2/nf-foo', '1.0.0')
        def manager = new EmbeddedPluginManager(store)

        when:
        manager.loadPlugins()

        then:
        manager.getPlugins()*.pluginId.sort() == ['nf-bar', 'org1/nf-baz', 'org1/nf-foo', 'org2/nf-foo']

        cleanup:
        store?.deleteDir()
    }
}
