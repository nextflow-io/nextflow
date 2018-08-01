/*
 * Copyright (c) 2018, University of Tübingen, Quantitative Biology Center (QBiC).
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
