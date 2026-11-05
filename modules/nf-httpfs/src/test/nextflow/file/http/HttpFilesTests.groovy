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

import java.nio.ByteBuffer
import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Paths

import com.github.tomakehurst.wiremock.junit.WireMockRule
import com.sun.net.httpserver.BasicAuthenticator
import com.sun.net.httpserver.Headers
import com.sun.net.httpserver.HttpExchange
import com.sun.net.httpserver.HttpHandler
import com.sun.net.httpserver.HttpServer
import groovy.transform.CompileStatic
import org.junit.Rule
import spock.lang.IgnoreIf
import spock.lang.Specification

import static com.github.tomakehurst.wiremock.client.WireMock.*

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HttpFilesTests extends Specification {

    @Rule
    WireMockRule wireMockRule = new WireMockRule(0)

    static final String BODY = 'line 1\nline 2\nline 3\n'

    static final String LAST_MODIFIED = 'Fri, 04 Nov 2016 21:50:34 GMT'

    private String stubFile(String name) {
        wireMockRule.stubFor(get(urlEqualTo(name))
            .willReturn(aResponse()
                .withStatus(200)
                .withBody(BODY)
                .withHeader('Content-Type', 'text/plain')
                .withHeader('Content-Length', BODY.size().toString())
                .withHeader('Last-Modified', LAST_MODIFIED)))
        return "http://localhost:${wireMockRule.port()}${name}"
    }

    private String stubMissing(String name) {
        wireMockRule.stubFor(get(urlEqualTo(name)).willReturn(aResponse().withStatus(404)))
        return "http://localhost:${wireMockRule.port()}${name}"
    }

    def 'should read a http file' () {
        given:
        def path = Paths.get(new URI(stubFile('/index.txt')))

        expect:
        path instanceof XPath
        path.fileSystem.provider() instanceof HttpFileSystemProvider
        and:
        path.text == BODY
        Files.readAllLines(path, Charset.forName('UTF-8')) == ['line 1', 'line 2', 'line 3']
        new String(Files.readAllBytes(path), Charset.forName('UTF-8')) == BODY
    }

    def 'should check file properties' () {
        given:
        def url = stubFile('/index.txt')
        def path1 = Paths.get(new URI(url))
        def path2 = Paths.get(new URI(stubMissing('/missing.txt')))
        def path3 = Paths.get(new URI(url))

        expect:
        Files.exists(path1)
        Files.size(path1) == BODY.size()
        Files.getLastModifiedTime(path1).toString() == '2016-11-04T21:50:34Z'
        !Files.isDirectory(path1)
        Files.isReadable(path1)
        !Files.isExecutable(path1)
        !Files.isWritable(path1)
        !Files.isHidden(path1)
        Files.isRegularFile(path1)
        !Files.isSymbolicLink(path1)
        Files.isSameFile(path1, path3)
        !Files.isSameFile(path1, path2)
        !Files.exists(path2)
    }

    def 'should copy a file' () {
        given:
        def source = Paths.get(new URI(stubFile('/index.txt')))
        def target = Files.createTempDirectory('nf').resolve('test.txt')

        when:
        Files.copy(source, target)
        then:
        target.text == BODY

        cleanup:
        target?.parent?.deleteDir()
    }

    def 'should re-try read file' () {
        given:
        def RESP = 'Hello world'
        and:
        // launch web server
        def attempt = 0
        HttpServer server = HttpServer.create(new InetSocketAddress(9900), 0);
        server.createContext("/", new BasicHandler(RESP, { ++attempt < 3 ? 404 : 200 } ));

        server.start()

        when:
        def path = Paths.get(new URI('http://admin:Secret1@localhost:9900/foo/bar'))
        then:
        path.text == 'Hello world'
        and:
        attempt == 3

        cleanup:
        server?.stop(0)
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    def 'should read FTP file' () {
        when:
        def lines = Paths.get(new URI('ftp://ftp.ncbi.nlm.nih.gov/robots.txt')).text.readLines()
        then:
        lines[0] == 'User-agent: *'
        lines[1] == 'Disallow: /'
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    def 'should read HTTPS file' () {
        given:
        def path = Paths.get(new URI('https://www.nextflow.io/index.html'))

        expect:
        path instanceof XPath
        path.fileSystem.provider() instanceof HttpsFileSystemProvider
        Files.exists(path)
        Files.readAllBytes(path).length > 0
    }

    def 'should read with a newByteChannel' () {
        given:
        def localhost = "http://localhost:${wireMockRule.port()}"
        wireMockRule.stubFor(get(urlEqualTo("/index.txt"))
            .willReturn(aResponse()
                .withStatus(200)
                .withBody("01234567890123456789012345678901234567890123456789")
                .withHeader("Content-Type", "text/html")
                .withHeader("Content-Length", "50")
                .withHeader("Last-Modified", "Fri, 04 Nov 2016 21:50:34 GMT")))

        when:
        def path = Paths.get(new URI("${localhost}/index.txt"))
        def sbc = Files.newByteChannel(path)
        def buffer = new ByteArrayOutputStream()
        ByteBuffer bf = ByteBuffer.allocate(15)
        while((sbc.read(bf))>0) {
            bf.flip();
            buffer.write(bf.array(), 0, bf.limit())
            bf.clear();
        }

        then:
        buffer.toString() == path.text
        buffer.toString() == '01234567890123456789012345678901234567890123456789'
    }

    def 'should use basic auth' () {
        given:
        def RESP = 'Hello world'
        and:
        // launch web server
        int port = new ServerSocket(0).withCloseable { it.localPort }
        HttpServer server = HttpServer.create(new InetSocketAddress(port), 0);
        def hc1 = server.createContext("/", new BasicHandler(RESP, 200));

        hc1.setAuthenticator(new BasicAuthenticator("get") {
            @Override
            boolean checkCredentials(String user, String pwd) {
                return user.equals("admin") && pwd.equals("Secret1");
            } });

        server.start()

        when:
        def path = Paths.get(new URI("http://admin:Secret1@localhost:${port}/foo/bar"))
        then:
        path.text == 'Hello world'

        cleanup:
        server?.stop(0)
    }

    @CompileStatic
    static class BasicHandler implements HttpHandler {

        String body

        Closure<Integer> respCode

        Map<String,String> allHeaders

        BasicHandler(String s, int code, Map<String,String> headers=null) {
            this.body=s
            this.respCode = { return code }
            this.allHeaders = headers ?: Collections.<String,String>emptyMap()
        }

        BasicHandler(String s, Closure<Integer> code, Map<String,String> headers=null) {
            this.body=s
            this.respCode = code
            this.allHeaders = headers ?: Collections.<String,String>emptyMap()
        }

        @Override
        void handle(HttpExchange request) throws IOException {

            Headers header = request.getResponseHeaders()

            for( Map.Entry<String,String> entry : allHeaders.entrySet() ) {
                header.set(entry.key, entry.value)
            }

            header.set("Content-Type", "text/plain")
            request.sendResponseHeaders(respCode.call(), body.size())

            OutputStream os = request.getResponseBody();
            os.write(body.bytes);
            os.close();
        }
    }
}
