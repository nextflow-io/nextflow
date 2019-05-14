/*
 * Copyright 2019, Genome Research Limited
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

package nextflow.wr.executor
import javax.net.ssl.HttpsURLConnection

import spock.lang.Specification
/**
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * based on K8sClientTest by Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WrRestApiTest extends Specification {

    // *** write tests

    // def 'should create a request' () {

    //     given:
    //     final TOKEN = '8d09d0ds'
    //     final client = Spy(WrRestApi)

    //     def HTTPS_CONN = Mock(HttpsURLConnection)
    //     def HTTP_CONN = Mock(HttpURLConnection)

    //     when:
    //     client.config.server = 'host.com:443'
    //     client.config.token = TOKEN
    //     def resp = client.makeRequest('GET', '/foo/bar')
    //     then:
    //     1 * client.createConnection0("https://host.com:443/foo/bar") >> HTTPS_CONN
    //     1 * client.setupHttpsConn(HTTPS_CONN) >> null
    //     1 * HTTPS_CONN.setRequestMethod('GET') >> null
    //     1 * HTTPS_CONN.setRequestProperty("Authorization", "Bearer $TOKEN")
    //     1 * HTTPS_CONN.setRequestProperty("Content-Type", "application/json")
    //     1 * HTTPS_CONN.getResponseCode() >> 200
    //     1 * HTTPS_CONN.getInputStream() >> { new ByteArrayInputStream('{"field_x":"hello"}'.bytes) }
    //     resp instanceof WrRestApi
    //     resp.text == '{"field_x":"hello"}'

    //     when:
    //     client.config.server = 'http://my-server.com'
    //     client.config.token = TOKEN
    //     client.makeRequest('POST', '/foo/bar')
    //     then:
    //     1 * client.createConnection0("http://my-server.com/foo/bar") >> HTTP_CONN
    //     0 * client.setupHttpsConn(_) >> null
    //     1 * HTTP_CONN.setRequestMethod('POST') >> null
    //     1 * HTTP_CONN.getResponseCode() >> 401
    //     1 * HTTP_CONN.getErrorStream() >> { new ByteArrayInputStream('{"field_x":"oops.."}'.bytes) }

    // }

    // def 'should make a get request' () {

    //     given:
    //     def client = Spy(WrRestApi)
    //     when:
    //     client.get('/foo/bar')
    //     then:
    //     1 * client.makeRequest('GET', '/foo/bar') >> null
    // }

    // def 'should make a post request' () {

    //     given:
    //     def client = Spy(WrRestApi)
    //     when:
    //     client.post('/foo/bar', '{ the: body }')
    //     then:
    //     1 * client.makeRequest('POST', '/foo/bar', '{ the: body }') >> null
    // }

    // def 'should make a delete request' () {

    //     given:
    //     def client = Spy(WrRestApi)
    //     when:
    //     client.delete('/foo/bar', '{ the: body }')
    //     then:
    //     1 * client.makeRequest('DELETE', '/foo/bar', '{ the: body }') >> null
    // }

}
