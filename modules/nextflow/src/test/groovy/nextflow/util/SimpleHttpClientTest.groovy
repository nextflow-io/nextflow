/*
 * Copyright 2018, University of Tübingen, Quantitative Biology Center (QBiC)
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.util

import com.sun.net.httpserver.HttpExchange
import com.sun.net.httpserver.HttpHandler
import com.sun.net.httpserver.HttpServer
import spock.lang.Specification
import spock.lang.Timeout

class SimpleHttpClientTest extends Specification{


    def 'should raise a ConnectException' () {

        given:
        def url = "http://localhost:3100"
        def httpClient = new SimpleHttpClient()

        when:
        httpClient.sendHttpMessage(url, '{"test_id": 2}')

        then:
        thrown(ConnectException)
        httpClient.getResponse() == null
        httpClient.getResponseCode() == -1

    }

    def 'should return a HTTP 404 error' () {

        given:
        def dummyUrl = "http://foo.bar"
        def con = Mock(HttpURLConnection)
        def httpClient = Spy(SimpleHttpClient)
        con.getOutputStream() >> new OutputStream() {
            @Override
            void write(int i) throws IOException { }
        }
        con.getInputStream() >> new InputStream() {
            @Override
            int read() throws IOException { return -1 }
        }
        con.getResponseCode() >> 404

        when:
        httpClient.sendHttpMessage(dummyUrl, '{"test_id": 2}')

        then:
        1 * httpClient.getHttpConnection(dummyUrl) >> con 
        httpClient.getResponseCode() == 404

    }


    def 'raise exception when host not available'() {

        given:
        def dummyUrl = "http://foo.bar"
        def httpClient = new SimpleHttpClient()

        when:
        httpClient.sendHttpMessage(dummyUrl, '{"test_id": "2"}')

        then:
        thrown(UnknownHostException)
    }

    def 'check isJson() method'() {

        when:
        def httpClient = new SimpleHttpClient()
        def noJson = httpClient.isJson("noJson")
        def json = httpClient.isJson('{"test": 2}')

        then:
        noExceptionThrown()
        !noJson
        json

    }

    def 'check for missing URL'() {
        given:
        def url = ''
        def httpClient = new SimpleHttpClient()

        when:
        httpClient.sendHttpMessage(url, '{"foo":, "bar"}')

        then:
        thrown(IllegalStateException)
    }

    @Timeout(1)
    def 'should make http request' () {

        given:
        def PAYLOAD = '{"hello":"world!"}'
        def RESULT = '{"status":"OK"}'
        def TOKEN = 'my:secret'
        def PORT = 9900
        def ENDPOINT = "http://localhost:$PORT/foo"
        def AGENT = 'NEXTFLOW/1.0'

        def expectedPayload
        def expectedContentType
        def expectedAuthToken
        def expectedUserAgent

        def server = HttpServer.create(new InetSocketAddress(PORT), 0)
        server.createContext("/", new HttpHandler() {
            @Override
            void handle(HttpExchange exchange) throws IOException {
                expectedPayload = exchange.requestBody.text
                expectedContentType = exchange.requestHeaders.get('Content-Type').get(0)
                expectedAuthToken = exchange.requestHeaders.get('Authorization').get(0)
                expectedUserAgent = exchange.requestHeaders.get('User-Agent').get(0)

                exchange.sendResponseHeaders(200, RESULT.getBytes().length)
                OutputStream os = exchange.getResponseBody()
                os.write(RESULT.getBytes());
                os.close()
            }
        });
        server.start()

        when:
        def client = new SimpleHttpClient()
        client.setAuthToken(TOKEN)
        client.setUserAgent(AGENT)
        client.sendHttpMessage( ENDPOINT, PAYLOAD )
        then:
        client.responseCode == 200
        client.response == RESULT
        and:
        expectedPayload == PAYLOAD
        expectedContentType == 'application/json'
        expectedAuthToken == "Basic ${TOKEN.bytes.encodeBase64()}"
        expectedUserAgent == AGENT

        cleanup:
        server?.stop(0)

    }


}
