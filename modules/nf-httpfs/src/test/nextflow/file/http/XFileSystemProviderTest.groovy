/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.file.http

import java.nio.file.Files
import java.nio.file.Path

import com.github.tomakehurst.wiremock.junit.WireMockRule
import com.github.tomjankes.wiremock.WireMockGroovy
import nextflow.SysEnv
import org.junit.Rule
import spock.lang.Specification
import spock.lang.Unroll
/**
 * Created by emilio on 08/11/16.
 */
class XFileSystemProviderTest extends Specification {

    def 'should create with default config settings' () {
        when:
        def fs = new HttpFileSystemProvider()
        then:
        fs.retryCodes() == HttpFileSystemProvider.DEFAULT_RETRY_CODES.tokenize(',').collect( it -> it as int )
        fs.backOffDelay() == HttpFileSystemProvider.DEFAULT_BACK_OFF_DELAY
        fs.backOffBase() == HttpFileSystemProvider.DEFAULT_BACK_OFF_BASE
        fs.maxAttempts() == HttpFileSystemProvider.DEFAULT_MAX_ATTEMPTS

    }

    def 'should create with custom config settings' () {
        given:
        SysEnv.push([NXF_HTTPFS_MAX_ATTEMPTS: '10',
                     NXF_HTTPFS_BACKOFF_BASE: '300',
                     NXF_HTTPFS_DELAY       : '400',
                     NXF_HTTPFS_RETRY_CODES : '1,2,3'])

        when:
        def fs = new HttpFileSystemProvider()
        then:
        fs.retryCodes() == [1,2,3]
        fs.backOffDelay() == 400
        fs.backOffBase() == 300
        fs.maxAttempts() == 10

        cleanup:
        SysEnv.pop()
    }

    def "should return input stream"() {
        given:
        def fsp = new HttpFileSystemProvider()
        def path = fsp.getPath(new URI('http://www.google.com/index.html'))
        when:
        def stream = fsp.newInputStream(path)
        then:
        stream.text.startsWith("<!doctype html>")
    }

    def "should return input stream from path"() {
        given:
        def fsp = Spy(new HttpFileSystemProvider())
        def path = fsp.getPath(new URI('http://host.com/index.html?query=123'))
        def connection = Mock(URLConnection)
        when:
        def stream = fsp.newInputStream(path)
        then:
        fsp.toConnection(path) >> { Path it ->
            assert it instanceof XPath
            assert it.toUri() == new URI('http://host.com/index.html?query=123')
            return connection
        }
        and:
        connection.getInputStream() >> new ByteArrayInputStream('Hello world'.bytes)
        and:
        stream.text == 'Hello world'
    }

    def "should read file attributes from map"() {
        given:
        def fs = new HttpFileSystemProvider()
        def attrMap = ['Last-Modified': ['Fri, 04 Nov 2016 21:50:34 GMT'], 'Content-Length': ['21729']]

        when:
        def attrs = fs.readHttpAttributes(attrMap)
        then:
        attrs.lastModifiedTime().toString() == '2016-11-04T21:50:34Z'
        attrs.size() == 21729

        when:
        attrs = fs.readHttpAttributes([:])
        then:
        attrs.lastModifiedTime() == null
        attrs.size() == -1
    }

    def "should read file attributes with german lang"() {
        given:
        def defLocale = Locale.getDefault(Locale.Category.FORMAT)
        // set german as current language
        def GERMAN = new Locale.Builder().setLanguage("de").setRegion("DE").build()
        Locale.setDefault(Locale.Category.FORMAT, GERMAN)
        def fs = new HttpFileSystemProvider()
        def attrMap = ['Last-Modified': ['Fri, 04 Nov 2016 21:50:34 GMT'], 'Content-Length': ['21729']]

        when:
        def attrs = fs.readHttpAttributes(attrMap)
        then:
        attrs.lastModifiedTime().toString() == '2016-11-04T21:50:34Z'
        attrs.size() == 21729

        cleanup:
        Locale.setDefault(Locale.Category.FORMAT, defLocale)
    }


    def "should read file attributes from HttpPath"() {
        given:
        def fsp = new HttpFileSystemProvider()
        def path = (XPath) fsp.getPath(new URI('http://www.nextflow.io/index.html'))

        when:
        def attrs = fsp.readHttpAttributes(path)
        then:
        attrs.lastModifiedTime() == null
        attrs.size() > 0
    }

    def "should read file attributes from FtpPath"() {
        given:
        def fsp = new FtpFileSystemProvider()
        def path = (XPath) fsp.getPath(new URI('ftp://ftp.ebi.ac.uk/robots.txt'))

        when:
        def attrs = fsp.readHttpAttributes(path)
        then:
        attrs.lastModifiedTime() == null
        attrs.size() == -1
    }

    @Unroll
    def 'should get uri path'() {
        given:
        def provider = new HttpFileSystemProvider()

        when:
        def path = provider.getPath(new URI(PATH))
        then:
        path.toUri().toString() == EXPECTED

        where:
        PATH                             | EXPECTED
        'http://foo.com/this/that'       | 'http://foo.com/this/that'
        'http://FOO.com/this/that'       | 'http://foo.com/this/that'
        'http://MrXYZ@foo.com/this/that' | 'http://MrXYZ@foo.com/this/that'
        'http://MrXYZ@FOO.com/this/that' | 'http://MrXYZ@foo.com/this/that'
        'http://@FOO.com/this/that'      | 'http://@foo.com/this/that'
        'http://foo.com/this/that?foo=1' | 'http://foo.com/this/that?foo=1'
    }

    @Unroll
    def 'should encode user info'() {
        given:
        def provider = new HttpFileSystemProvider()
        expect:
        provider.auth(USER_INFO) == EXPECTED
        where:
        USER_INFO              | EXPECTED
        "foo:bar"              | "Basic ${'foo:bar'.bytes.encodeBase64()}"
        "x-oauth-bearer:12345" | "Bearer 12345"
    }

    @Rule
    WireMockRule wireMockRule = new WireMockRule(18080)

    def 'should follow a redirect when read a http file '() {
        given:
        def wireMock = new WireMockGroovy(18080)
        wireMock.stub {
            request {
                method "GET"
                url "/index.html"
            }
            response {
                status HTTP_CODE
                headers {
                    "Location" "http://localhost:18080${REDIRECT_TO}"
                }
            }
        }
        wireMock.stub {
            request {
                method "GET"
                url "/index2.html"
            }
            response {
                status HTTP_CODE
                headers {
                    "Location" "http://localhost:18080/redirected.html"
                }
            }
        }
        wireMock.stub {
            request {
                method "GET"
                url "/redirected.html"
            }
            response {
                status 200
                body """a
                 b
                 c
                 d
                 """
                headers {
                    "Content-Type" "text/html"
                    "Content-Length" "10"
                    "Last-Modified" "Fri, 04 Nov 2016 21:50:34 GMT"
                }
            }
        }
        and:
        def provider = new HttpFileSystemProvider()
        when:
        def path = provider.getPath(new URI('http://localhost:18080/index.html'))
        then:
        path
        Files.size(path) == EXPECTED

        where:
        HTTP_CODE | REDIRECT_TO        | EXPECTED
        300       | "/redirected.html" | 10
        300       | "/index2.html"     | 10

        301       | "/redirected.html" | 10
        301       | "/index2.html"     | 10

        302       | "/redirected.html" | 10
        302       | "/index2.html"     | 10

        307       | "/redirected.html" | 10
        307       | "/index2.html"     | 10

        308       | "/redirected.html" | 10
        308       | "/index2.html"     | 10
        //infinite redirect to himself
        308       | "/index.html"      | -1
    }
}
