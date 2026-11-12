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
import java.nio.file.Path

import spock.lang.Specification

class LocalPluginManagerTest extends Specification {

    static Path createPlugin(Path store, String id, String ver) {
        final dir = store.resolve("$id-$ver")
        dir.mkdirs()
        dir.resolve('MANIFEST.MF').text = """\
            Manifest-Version: 1.0
            Plugin-Class: ${PluginUpdaterTest.FooPlugin.name}
            Plugin-Id: $id
            Plugin-Version: $ver
            """.stripIndent()
        return dir
    }

    def 'should create a unique link for plugins with the same name in different namespaces' () {
        given:
        def store = Files.createTempDirectory('test')
        def manager = new LocalPluginManager(store)

        when:
        manager.loadPluginFromPath(createPlugin(store, 'org1/nf-foo', '1.0.0'))
        manager.loadPluginFromPath(createPlugin(store, 'org2/nf-foo', '1.0.0'))
        manager.loadPluginFromPath(createPlugin(store, 'nf-bar', '1.0.0'))

        then:
        manager.localRoot.resolve('org1~nf-foo-1.0.0').isLink()
        manager.localRoot.resolve('org2~nf-foo-1.0.0').isLink()
        manager.localRoot.resolve('nf-bar-1.0.0').isLink()
        manager.getPlugin('org1/nf-foo') != null
        manager.getPlugin('org2/nf-foo') != null

        cleanup:
        store?.deleteDir()
    }
}
