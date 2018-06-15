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

import spock.lang.Shared
import spock.lang.Specification

class SimpleHttpClientTest extends Specification{

    @Shared SimpleHttpClient httpClient

    @Shared String dummyUrl

    def setupSpec(){

        httpClient = new SimpleHttpClient()

        dummyUrl = "http://foo.bar"

    }

    def 'should set a http url successfully' () {

        given:
        def url = "http://localhost"

        when:
        httpClient.setUrl(url)

        then:
        noExceptionThrown()
        assert httpClient.getUrl() == url

    }

    def 'should raise a ConnectException' () {

        given:
        def url = "http://localhost"

        when:
        httpClient.setUrl(url)
        httpClient.sendHttpMessage('{"test_id": 2}')

        then:
        thrown(ConnectException)
        assert httpClient.getResponse() == ''
        assert httpClient.getResponseCode() == -1

    }

    def 'should return a HTTP 404 error' () {

        given:
        def con = Mock(HttpURLConnection)
        def httpClient0 = Spy(SimpleHttpClient)
        httpClient0.setUpConnection(dummyUrl) >> con
        httpClient0.setUrl(dummyUrl)
        con.getOutputStream() >> new OutputStream() {
            @Override
            void write(int i) throws IOException {

            }
        }
        con.getInputStream() >> new InputStream() {
            @Override
            int read() throws IOException {
                return -1
            }
        }
        con.getResponseCode() >> 404

        when:
        httpClient0.sendHttpMessage('{"test_id": 2}')

        then:
        noExceptionThrown()
        assert httpClient0.getResponseCode() == 404 : "Response code was ${httpClient.getResponseCode()}"

    }

    def 'should not send message when not a JSON-like string' (){

        given:
        def httpClient0 = Spy(SimpleHttpClient)
        httpClient0.setUrl(dummyUrl)

        when:
        httpClient0.sendHttpMessage("Foobar")

        then:
        thrown(IllegalArgumentException)
        assert httpClient0.getResponseCode() == -1
        1 * httpClient0.isJson("Foobar")
        0 * httpClient0.setUpConnection(dummyUrl)

    }

    def 'raise exception when host not available'() {

        given:
        httpClient.setUrl(dummyUrl)

        when:
        httpClient.sendHttpMessage('{"test_id": 2}')

        then:
        thrown(UnknownHostException)
    }

    def 'check isJson() method'() {

        when:
        def noJson = httpClient.isJson("noJson")
        def json = httpClient.isJson('{"test": 2}')

        then:
        noExceptionThrown()
        assert !noJson
        assert json

    }

    def 'ckeck for missing URL'() {

        when:
        httpClient.setUrl("")
        httpClient.sendHttpMessage()


        then:
        thrown(IllegalStateException)


    }


}
