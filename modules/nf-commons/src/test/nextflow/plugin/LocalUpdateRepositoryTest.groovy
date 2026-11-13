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

class LocalUpdateRepositoryTest extends Specification {

    def 'should find plugins in namespace folders' () {
        given:
        def store = Files.createTempDirectory('test')
        createPlugin(store, 'nf-bar', '1.0.0')
        createPlugin(store, 'org1/nf-foo', '1.0.0')
        createPlugin(store, 'org1/nf-foo', '2.0.0')
        createPlugin(store, 'org2/nf-foo', '1.0.0')
        def repo = new LocalUpdateRepository('local', store)

        expect:
        repo.getPlugins().keySet() == ['nf-bar', 'org1/nf-foo', 'org2/nf-foo'] as Set
        repo.getPlugin('org1/nf-foo').releases*.version.sort() == ['1.0.0', '2.0.0']

        cleanup:
        store?.deleteDir()
    }
}
