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

package nextflow.util


import spock.lang.Specification 
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProxyConfigTest extends Specification {


    def 'should parse proxy env variables'( ) {

        expect:
        ProxyConfig.parse(null) == null

        ProxyConfig.parse('http://domain') == new ProxyConfig(protocol: 'http', host: 'domain')
        ProxyConfig.parse('http://domain:333') == new ProxyConfig(protocol: 'http',host: 'domain', port: '333')
        ProxyConfig.parse('http://10.20.30.40') == new ProxyConfig(protocol: 'http',host: '10.20.30.40')
        ProxyConfig.parse('http://10.20.30.40:333') == new ProxyConfig(protocol: 'http',host: '10.20.30.40', port: '333')
        ProxyConfig.parse('http://10.20.30.40:333/some/path') == new ProxyConfig(protocol: 'http',host: '10.20.30.40', port: '333')

        ProxyConfig.parse('http://user:pass@domain') == new ProxyConfig(protocol: 'http',host: 'domain', username: 'user', password: 'pass')
        ProxyConfig.parse('http://user:pass@domain:333') == new ProxyConfig(protocol: 'http',host: 'domain', port: '333', username: 'user', password: 'pass')
        ProxyConfig.parse('http://user:pass@10.20.30.40') == new ProxyConfig(protocol: 'http',host: '10.20.30.40', username: 'user', password: 'pass')
        ProxyConfig.parse('http://user:pass@10.20.30.40:333') == new ProxyConfig(protocol: 'http',host: '10.20.30.40', port: '333', username: 'user', password: 'pass')
        ProxyConfig.parse('http://user:pass@10.20.30.40:333/some/path') == new ProxyConfig(protocol: 'http',host: '10.20.30.40', port: '333', username: 'user', password: 'pass')

        ProxyConfig.parse('foo') == new ProxyConfig(host: 'foo')
        ProxyConfig.parse('foo:123') == new ProxyConfig(host: 'foo', port: '123')

    }
    
}
