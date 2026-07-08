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

import org.pf4j.update.PluginInfo
import spock.lang.Specification
import test.EnvHelper
import test.MockAuthProxyServer

/**
 * Verify that HTTP clients created via a plain {@code HxClient.newBuilder()}
 * (e.g. the plugin registry client) automatically resolve the forward proxy
 * from the environment, including proxy credentials.
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
class HttpPluginRepositoryProxyTest extends Specification {

    def 'should fetch plugin metadata via an authenticated proxy'() {
        given: 'an authenticating forward proxy answering with plugin metadata'
        def proxy = new MockAuthProxyServer('foo', 'secret').start()
        proxy.responseContentType = 'application/json'
        proxy.responseBody = '''{
              "plugins": [
                {
                  "id": "nf-fake",
                  "releases": [
                    {
                      "version": "0.0.1",
                      "url": "http://example.com/fake-plugin/0.0.1/plugin.zip",
                      "date": "2025-03-21T10:40:36Z",
                      "sha512sum": "cbc4a7f0bc10c955ff10b85da85f4df8303e4174f1c44922c52a274a59afd786ee3bc685f332d1e0a7497d0b28b0cfc9939d052d6fa992607f43c870825f6caf",
                      "requires": ">=24.10.0"
                    }
                  ]
                }
              ]
            }'''
        when: 'the registry host is only reachable via the proxy'
        PluginInfo info = null
        EnvHelper.withEnv(proxy.proxyEnv()) {
            final repo = new HttpPluginRepository('test-repo', new URI('http://plugins.registry.internal/'))
            info = repo.getPlugin('nf-fake')
        }

        then:
        info.id == 'nf-fake'
        info.releases.size() == 1
        info.releases[0].version == '0.0.1'
        and: 'the request went through the proxy after the 407 challenge'
        proxy.unauthorizedCount.get() >= 1
        proxy.proxiedCount.get() >= 1

        cleanup:
        proxy?.close()
    }
}
