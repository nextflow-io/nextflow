/*
 * Copyright 2018, University of Tübingen, Quantitative Biology Center (QBiC)
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import spock.lang.Specification

class SimpleHttpClientTest extends Specification{


    def 'should set a http url successfully' () {

        given:
        def url = "http://localhost"
        def httpClient = new SimpleHttpClient()
        when:
        httpClient.setUrl(url)

        then:
        noExceptionThrown()
        httpClient.getUrl() == url
    }


    def 'should raise a ConnectException' () {

        given:
        def url = "http://localhost:3100"
        def httpClient = new SimpleHttpClient()

        when:
        httpClient.setUrl(url)
        httpClient.sendHttpMessage('{"test_id": 2}')

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
        httpClient.setUpConnection(dummyUrl) >> con
        httpClient.setUrl(dummyUrl)
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
        httpClient.sendHttpMessage('{"test_id": 2}')

        then:
        noExceptionThrown()
        httpClient.getResponseCode() == 404

    }

    def 'should not send message when not a JSON-like string' (){

        given:
        def dummyUrl = "http://foo.bar"
        def httpClient0 = Spy(SimpleHttpClient)
        httpClient0.setUrl(dummyUrl)

        when:
        httpClient0.sendHttpMessage("Foobar")

        then:
        thrown(IllegalArgumentException)
        httpClient0.getResponseCode() == -1
        1 * httpClient0.isJson("Foobar")
        0 * httpClient0.setUpConnection(dummyUrl)

    }

    def 'raise exception when host not available'() {

        given:
        def dummyUrl = "http://foo.bar"
        def httpClient = new SimpleHttpClient()
        httpClient.setUrl(dummyUrl)

        when:
        httpClient.sendHttpMessage('{"test_id": "2"}')

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
        def httpClient = new SimpleHttpClient()

        when:
        httpClient.setUrl("")
        httpClient.sendHttpMessage()

        then:
        thrown(IllegalStateException)
    }

}
