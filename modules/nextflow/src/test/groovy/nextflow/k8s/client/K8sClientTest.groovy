/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.k8s.client

import javax.net.ssl.HttpsURLConnection

import nextflow.k8s.client.res.PodProvider
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sClientTest extends Specification {

    def 'should create a request' () {

        given:
        final TOKEN = '8d09d0ds'
        final client = Spy(K8sRequestDelegate)

        def HTTPS_CONN = Mock(HttpsURLConnection)
        def HTTP_CONN = Mock(HttpURLConnection)

        when:
        client.config.server = 'host.com:443'
        client.config.token = TOKEN
        def resp = client.makeRequest('GET', '/foo/bar')
        then:
        1 * client.createConnection0("https://host.com:443/foo/bar") >> HTTPS_CONN
        1 * client.setupHttpsConn(HTTPS_CONN) >> null
        1 * HTTPS_CONN.setRequestMethod('GET') >> null
        1 * HTTPS_CONN.setRequestProperty("Authorization", "Bearer $TOKEN")
        1 * HTTPS_CONN.setRequestProperty("Content-Type", "application/json")
        1 * HTTPS_CONN.getResponseCode() >> 200
        1 * HTTPS_CONN.getInputStream() >> { new ByteArrayInputStream('{"field_x":"hello"}'.bytes) }
        resp instanceof K8sResponseApi
        resp.text == '{"field_x":"hello"}'

        when:
        client.config.server = 'http://my-server.com'
        client.config.token = TOKEN
        client.makeRequest('POST', '/foo/bar')
        then:
        1 * client.createConnection0("http://my-server.com/foo/bar") >> HTTP_CONN
        0 * client.setupHttpsConn(_) >> null
        1 * HTTP_CONN.setRequestMethod('POST') >> null
        1 * HTTP_CONN.getResponseCode() >> 401
        1 * HTTP_CONN.getErrorStream() >> { new ByteArrayInputStream('{"field_x":"oops.."}'.bytes) }
        def e = thrown(K8sResponseException)
        e.response.field_x == 'oops..'

    }

    def 'should delegate resource requests to provider' () {
        given:
        def result
        def provider = Spy(PodProvider)
        def client = Spy(new K8sClient(new ClientConfig(), ['Pod': provider]))
        def RESP = [response: 'hello']
        and:
        def spec = [foo: 'bar']
        def stream = Mock(InputStream)

        when:
        result = client.create('Pod', spec)
        then:
        1 * provider.create(spec, null) >> RESP
        result.response == 'hello'

        when:
        result = client.delete('Pod', 'foo')
        then:
        1 * provider.delete('foo') >> RESP
        result.response == 'hello'

        when:
        result = client.list('Pod')
        then:
        1 * provider.list(false) >> RESP
        result.response == 'hello'

        when:
        result = client.list('Pod', true)
        then:
        1 * provider.list(true) >> RESP
        result.response == 'hello'

        when:
        result = client.fetchLog('Pod', 'foo')
        then:
        1 * provider.fetchLog([:], 'foo') >> stream
        result == stream

        when:
        result = client.getState('Pod', 'foo')
        then:
        1 * provider.getState('foo') >> [status: 'COMPLETE']
        result.status == 'COMPLETE'

        when:
        result = client.podHostname('foo')
        then:
        1 * provider.getHostname('foo') >> 'hostname'
        result == 'hostname'
    }

    def 'should list secrets' () {
        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"hello"}'

        when:
        result = client.secretList()
        then:
        1 * client.makeRequest('GET',"/api/v1/namespaces/default/secrets") >> RESP
        result.response == 'hello'

    }

    def 'should describe secret' () {
        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"hello"}'

        when:
        result = client.secretDescribe('my-secret')
        then:
        1 * client.makeRequest('GET',"/api/v1/namespaces/default/secrets/my-secret") >> RESP
        result.response == 'hello'

        when:
        client.config.namespace = 'paperino'
        result = client.secretDescribe('pluto')
        then:
        1 * client.makeRequest('GET','/api/v1/namespaces/paperino/secrets/pluto') >> RESP
        result.response == 'hello'
    }

    def 'should create config' () {

        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"done"}'
        def CONFIG = [foo: 'hello', bar:'world']
        def JSON1 = '{"apiVersion":"v1","kind":"ConfigMap","metadata":{"name":"foo","namespace":"default"},"data":{"foo":"hello","bar":"world"}}'
        def JSON2 = '{"apiVersion":"v1","kind":"ConfigMap","metadata":{"name":"foo","namespace":"bar"},"data":{"foo":"hello","bar":"world"}}'

        when:
        result = client.configCreate('foo' , CONFIG)
        then:
        1 * client.makeRequest('POST', "/api/v1/namespaces/default/configmaps", JSON1) >> RESP
        result.response == 'done'

        when:
        client.config.namespace = 'bar'
        result = client.configCreate('foo' , CONFIG)
        then:
        1 * client.makeRequest('POST', "/api/v1/namespaces/bar/configmaps", JSON2) >> RESP
        result.response == 'done'

    }

    def 'should delete config' () {

        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"done"}'

        when:
        result = client.configDelete('foo')
        then:
        1 * client.makeRequest('DELETE',"/api/v1/namespaces/default/configmaps/foo") >> RESP
        result.response == 'done'

        when:
        client.config.namespace = 'ns-1'
        result = client.configDelete('foo')
        then:
        1 * client.makeRequest('DELETE',"/api/v1/namespaces/ns-1/configmaps/foo") >> RESP
        result.response == 'done'

    }

    def 'should delete all configs' () {

        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"done"}'

        when:
        result = client.configDeleteAll()
        then:
        1 * client.makeRequest('DELETE',"/api/v1/namespaces/default/configmaps") >> RESP
        result.response == 'done'

        when:
        client.config.namespace = 'ns-1'
        result = client.configDeleteAll()
        then:
        1 * client.makeRequest('DELETE',"/api/v1/namespaces/ns-1/configmaps") >> RESP
        result.response == 'done'

    }

}
