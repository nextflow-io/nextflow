/*
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

import java.nio.ByteBuffer
import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Paths

import com.github.tomakehurst.wiremock.junit.WireMockRule
import com.github.tomjankes.wiremock.WireMockGroovy
import org.junit.Rule
import spock.lang.IgnoreIf
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HttpFilesTests extends Specification {

    @Rule
    WireMockRule wireMockRule = new WireMockRule(18080)

    def 'should read http file from WireMock' () {

        given:
        def wireMock = new WireMockGroovy(18080)
        wireMock.stub {
            request {
                method "GET"
                url "/index.html"
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

        when:
        def path = Paths.get(new URI('http://localhost:18080/index.html'))
        then:
        Files.size(path) == 10
        Files.getLastModifiedTime(path).toString() == "2016-11-04T21:50:34Z"

        when:
        def path2 = Paths.get(new URI('http://localhost:18080/missing.html'))
        then:
        !Files.exists(path2)
    }


    def 'read a http file ' () {
        given:
        def uri = new URI('http://www.nextflow.io/index.html')
        when:
        def path = Paths.get(uri)
        then:
        path instanceof XPath

        when:
        def lines = Files.readAllLines(path, Charset.forName('UTF-8'))
        then:
        lines.size()>0
        lines[0] == '<html>'

    }

    def 'should check file properties' () {

        when:
        def path1 = Paths.get(new URI('http://www.nextflow.io/index.html'))
        def path2 = Paths.get(new URI('http://www.google.com/unknown'))
        def path3 = Paths.get(new URI('http://www.nextflow.io/index.html'))

        then:
        Files.exists(path1)
        Files.size(path1) > 0
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

    @IgnoreIf({ env['TRAVIS'] })
    def 'should read FTP file' () {
        when:
        def lines = Paths.get(new URI('ftp://ftp.ebi.ac.uk/robots.txt')).text.readLines()
        then:
        lines[0] == 'User-agent: *'
        lines[1] == 'Disallow: /'
    }

    def 'should read HTTPS file' () {

        given:
        def uri = new URI('https://www.nextflow.io/index.html')
        when:
        def path = Paths.get(uri)
        then:
        path instanceof XPath

        when:
        def lines = Files.readAllLines(path, Charset.forName('UTF-8'))
        then:
        lines.size()>0
        lines[0] == '<!DOCTYPE html>'

    }

    def 'should copy a file' () {

        given:
        def uri = new URI('https://www.nextflow.io/index.html')
        def source = Paths.get(uri)
        def target = Files.createTempDirectory('nf').resolve('test.html')

        when:
        Files.copy(source, target)

        then:
        source.text == target.text

        cleanup:
        target?.parent?.deleteDir()

    }

    @IgnoreIf({ env['TRAVIS'] })
    def 'should read lines' () {
        given:
        def path = Paths.get(new URI('ftp://ftp.ebi.ac.uk/robots.txt'))

        when:
        def lines = Files.readAllLines(path, Charset.forName('UTF-8'))
        then:
        lines[0] == 'User-agent: *'
        lines[1] == 'Disallow: /'
    }

    @IgnoreIf({ env['TRAVIS'] })
    def 'should read all bytes' ( ) {
        given:
        def path = Paths.get(new URI('ftp://ftp.ebi.ac.uk/robots.txt'))

        when:
        def bytes = Files.readAllBytes(path)
        def lines = new String(bytes).readLines()
        then:
        lines[0] == 'User-agent: *'
        lines[1] == 'Disallow: /'

    }

    def 'should read with a newByteChannel' () {

        given:
        def wireMock = new WireMockGroovy(18080)
        wireMock.stub {
            request {
                method "GET"
                url "/index.txt"
            }
            response {
                status 200
                body "01234567890123456789012345678901234567890123456789"
                headers {
                    "Content-Type" "text/html"
                    "Content-Length" "50"
                    "Last-Modified" "Fri, 04 Nov 2016 21:50:34 GMT"
                }
            }
        }

        when:
        def path = Paths.get(new URI('http://localhost:18080/index.txt'))
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

}
