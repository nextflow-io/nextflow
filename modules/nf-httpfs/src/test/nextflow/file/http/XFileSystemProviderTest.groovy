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

package nextflow.file.http

import java.nio.file.FileAlreadyExistsException
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption

import com.github.tomakehurst.wiremock.junit.WireMockRule
import nextflow.file.HttpCopyOption
import org.junit.Rule
import spock.lang.IgnoreIf
import spock.lang.Specification
import spock.lang.Unroll

import static com.github.tomakehurst.wiremock.client.WireMock.*
/**
 * Created by emilio on 08/11/16.
 */
class XFileSystemProviderTest extends Specification {

    @IgnoreIf({System.getenv('NXF_SMOKE')})
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
        def DATA = 'Hello world'
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
        connection.getInputStream() >> new ByteArrayInputStream(DATA.bytes)
        connection.getContentLengthLong() >> DATA.size()
        and:
        stream.text == 'Hello world'
    }

    def "should read file attributes from map"() {
        given:
        def fs = new HttpFileSystemProvider()
        def attrMap = ['last-modified': ['Fri, 04 Nov 2016 21:50:34 GMT'], 'content-length': ['21729']]

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
        def attrMap = ['last-modified': ['Fri, 04 Nov 2016 21:50:34 GMT'], 'content-length': ['21729']]

        when:
        def attrs = fs.readHttpAttributes(attrMap)
        then:
        attrs.lastModifiedTime().toString() == '2016-11-04T21:50:34Z'
        attrs.size() == 21729

        cleanup:
        Locale.setDefault(Locale.Category.FORMAT, defLocale)
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    def "should read file attributes from HttpPath"() {
        given:
        def fsp = new HttpFileSystemProvider()
        def path = (XPath) fsp.getPath(new URI('http://www.nextflow.io/index.html'))

        when:
        def attrs = fsp.readHttpAttributes(path)
        then:
        attrs.lastModifiedTime()
        attrs.size() > 0
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
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
    WireMockRule wireMockRule = new WireMockRule(0)

    @Unroll
    def 'should follow a redirect when read a http file '() {
        given:
        def localhost = "http://localhost:${wireMockRule.port()}"
        wireMockRule.stubFor(get(urlEqualTo("/index.html"))
            .willReturn(aResponse()
                .withStatus(HTTP_CODE)
                .withHeader("Location", "${localhost}${REDIRECT_TO}")))

        wireMockRule.stubFor(get(urlEqualTo("/index2.html"))
            .willReturn(aResponse()
                .withStatus(HTTP_CODE)
                .withHeader("Location", "${localhost}/target.html")))

        wireMockRule.stubFor(get(urlEqualTo("/target.html"))
            .willReturn(aResponse()
                .withStatus(200)
                .withBody("""a
                 b
                 c
                 d
                 """)
                .withHeader("Content-Type", "text/html")
                .withHeader("Content-Length", "10")
                .withHeader("Last-Modified", "Fri, 04 Nov 2016 21:50:34 GMT")))

        and:
        def provider = new HttpFileSystemProvider()
        when:
        def path = provider.getPath(new URI("${localhost}/index.html"))
        then:
        path
        Files.size(path) == EXPECTED

        where:
        HTTP_CODE | REDIRECT_TO         | EXPECTED
        301       | "/target.html"      | 10
        301       | "/index2.html"      | 10

        302       | "/target.html"      | 10
        302       | "/index2.html"      | 10

        307       | "/target.html"      | 10
        307       | "/index2.html"      | 10

        308       | "/target.html"      | 10
        308       | "/index2.html"      | 10
        //infinite redirect to himself
        308       | "/index.html"       | -1
    }

    def 'should normalize location' () {
        given:
        def provider = Spy(XFileSystemProvider)

        expect:
        provider.absLocation(LOCATION, new URL(TARGET)) == EXPECTED

        where:
        LOCATION                | TARGET                    | EXPECTED
        'https://this/that'     | 'http://foo.com:123'      | 'https://this/that'
        '/'                     | 'http://foo.com:123'      | 'http://foo.com:123/'
        '/this/that'            | 'http://foo.com:123'      | 'http://foo.com:123/this/that'
        '/this/that'            | 'http://foo.com:123/abc'  | 'http://foo.com:123/this/that'
        'this/that'             | 'http://foo.com:123/abc'  | 'http://foo.com:123/this/that'

    }

    def 'should resume an http download using a byte range'() {
        given:
        def localhost = "http://localhost:${wireMockRule.port()}"
        def FULL = 'ABCDEFGHIJKLMNOPQRST'   // 20 bytes
        def PARTIAL = 'ABCDEFGHIJ'          // first 10 bytes
        def REMAINING = 'KLMNOPQRST'        // remaining 10 bytes

        wireMockRule.stubFor(get(urlEqualTo('/file.txt'))
            .withHeader('Range', equalTo('bytes=10-'))
            .willReturn(aResponse()
                .withStatus(206)
                .withHeader('Content-Range', 'bytes 10-19/20')
                .withHeader('Content-Length', '10')
                .withBody(REMAINING)))

        def provider = new HttpFileSystemProvider()
        def source = provider.getPath(new URI("${localhost}/file.txt"))
        def target = Files.createTempFile('nf-resume', '.txt')
        target.text = PARTIAL

        when:
        provider.download(source, target, HttpCopyOption.RESUME)

        then:
        target.text == FULL

        and:
        wireMockRule.verify(getRequestedFor(urlEqualTo('/file.txt')).withHeader('Range', equalTo('bytes=10-')))

        cleanup:
        target?.delete()
    }

    def 'should restart an http download when the server ignores the range'() {
        given:
        def localhost = "http://localhost:${wireMockRule.port()}"
        def FULL = 'ABCDEFGHIJKLMNOPQRST'
        def PARTIAL = 'ABCDEFGHIJ'

        wireMockRule.stubFor(get(urlEqualTo('/file.txt'))
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Length', '20')
                .withBody(FULL)))

        def provider = new HttpFileSystemProvider()
        def source = provider.getPath(new URI("${localhost}/file.txt"))
        def target = Files.createTempFile('nf-restart', '.txt')
        target.text = PARTIAL

        when:
        provider.download(source, target, HttpCopyOption.RESUME)

        then:
        target.text == FULL

        and:
        wireMockRule.verify(getRequestedFor(urlEqualTo('/file.txt')).withHeader('Range', equalTo('bytes=10-')))

        cleanup:
        target?.delete()
    }

    def 'should gate resume download to http and https'() {
        given:
        def http = new HttpFileSystemProvider()
        def https = new HttpsFileSystemProvider()
        def ftp = new FtpFileSystemProvider()

        expect:
        http.canDownload(null, null)
        https.canDownload(null, null)
        !ftp.canDownload(null, null)

        and:
        !http.canUpload(null, null)
        !https.canUpload(null, null)
        !ftp.canUpload(null, null)
    }

    def 'should restart when the server returns a mismatched content range'() {
        given:
        def localhost = "http://localhost:${wireMockRule.port()}"
        def FULL = 'ABCDEFGHIJKLMNOPQRST'   // 20 bytes
        def PARTIAL = 'ABCDEFGHIJ'          // first 10 bytes

        // ranged request returns a 206 that does not start at the requested offset
        wireMockRule.stubFor(get(urlEqualTo('/file.txt'))
            .atPriority(1)
            .withHeader('Range', equalTo('bytes=10-'))
            .willReturn(aResponse()
                .withStatus(206)
                .withHeader('Content-Range', 'bytes 0-19/20')
                .withHeader('Content-Length', '20')
                .withBody(FULL)))

        // fallback request returns the full body
        wireMockRule.stubFor(get(urlEqualTo('/file.txt'))
            .atPriority(2)
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Length', '20')
                .withBody(FULL)))

        def provider = new HttpFileSystemProvider()
        def source = provider.getPath(new URI("${localhost}/file.txt"))
        def target = Files.createTempFile('nf-mismatch', '.txt')
        target.text = PARTIAL

        when:
        provider.download(source, target, HttpCopyOption.RESUME)

        then:
        target.text == FULL

        and:
        wireMockRule.verify(getRequestedFor(urlEqualTo('/file.txt')).withHeader('Range', equalTo('bytes=10-')))

        cleanup:
        target?.delete()
    }

    def 'should restart when the server returns 416 range not satisfiable'() {
        given:
        def localhost = "http://localhost:${wireMockRule.port()}"
        def FULL = 'ABCDEFGHIJKLMNOPQRST'
        def PARTIAL = 'ABCDEFGHIJ'

        wireMockRule.stubFor(get(urlEqualTo('/file.txt'))
            .atPriority(1)
            .withHeader('Range', equalTo('bytes=10-'))
            .willReturn(aResponse()
                .withStatus(416)
                .withHeader('Content-Range', 'bytes */20')))

        wireMockRule.stubFor(get(urlEqualTo('/file.txt'))
            .atPriority(2)
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Length', '20')
                .withBody(FULL)))

        def provider = new HttpFileSystemProvider()
        def source = provider.getPath(new URI("${localhost}/file.txt"))
        def target = Files.createTempFile('nf-416', '.txt')
        target.text = PARTIAL

        when:
        provider.download(source, target, HttpCopyOption.RESUME)

        then:
        target.text == FULL

        and:
        wireMockRule.verify(getRequestedFor(urlEqualTo('/file.txt')).withHeader('Range', equalTo('bytes=10-')))

        cleanup:
        target?.delete()
    }

    def 'should download a file when the target does not exist'() {
        given:
        def localhost = "http://localhost:${wireMockRule.port()}"
        def FULL = 'ABCDEFGHIJKLMNOPQRST'

        wireMockRule.stubFor(get(urlEqualTo('/file.txt'))
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Length', '20')
                .withBody(FULL)))

        def provider = new HttpFileSystemProvider()
        def source = provider.getPath(new URI("${localhost}/file.txt"))
        def target = Files.createTempFile('nf-fresh', '.txt')
        Files.delete(target)

        when:
        provider.download(source, target)

        then:
        target.text == FULL

        cleanup:
        target?.delete()
    }

    def 'should replace an existing file when REPLACE_EXISTING is specified'() {
        given:
        def localhost = "http://localhost:${wireMockRule.port()}"
        def FULL = 'ABCDEFGHIJKLMNOPQRST'

        wireMockRule.stubFor(get(urlEqualTo('/file.txt'))
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Length', '20')
                .withBody(FULL)))

        def provider = new HttpFileSystemProvider()
        def source = provider.getPath(new URI("${localhost}/file.txt"))
        def target = Files.createTempFile('nf-replace', '.txt')
        target.text = 'stale content'

        when:
        provider.download(source, target, StandardCopyOption.REPLACE_EXISTING)

        then:
        target.text == FULL

        cleanup:
        target?.delete()
    }

    def 'should fail when the target already exists without REPLACE_EXISTING'() {
        given:
        def localhost = "http://localhost:${wireMockRule.port()}"
        def FULL = 'ABCDEFGHIJKLMNOPQRST'

        wireMockRule.stubFor(get(urlEqualTo('/file.txt'))
            .willReturn(aResponse()
                .withStatus(200)
                .withHeader('Content-Length', '20')
                .withBody(FULL)))

        def provider = new HttpFileSystemProvider()
        def source = provider.getPath(new URI("${localhost}/file.txt"))
        def target = Files.createTempFile('nf-exists', '.txt')
        target.text = 'stale content'

        when:
        provider.download(source, target)

        then:
        thrown(FileAlreadyExistsException)
        target.text == 'stale content'

        cleanup:
        target?.delete()
    }
}
