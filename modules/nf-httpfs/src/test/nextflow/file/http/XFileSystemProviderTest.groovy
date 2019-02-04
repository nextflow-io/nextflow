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
import spock.lang.Specification
/**
 * Created by emilio on 08/11/16.
 */
class XFileSystemProviderTest extends Specification {

    def "should return input stream"() {
        given:
        def fsp = new HttpFileSystemProvider()
        def path = fsp.getPath(new URI('http://www.google.com/index.html'))
        when:
        def stream = fsp.newInputStream(path)
        then:
        stream.text.startsWith("<!doctype html>")
    }

    def "should read file attributes from map"() {
        given:
        def fs = new HttpFileSystemProvider()
        def attrMap = ['Last-Modified': ['Fri, 04 Nov 2016 21:50:34 GMT'], 'Content-Length': ['21729'] ]

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
        def attrMap = ['Last-Modified': ['Fri, 04 Nov 2016 21:50:34 GMT'], 'Content-Length': ['21729'] ]

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
        def path = (XPath)fsp.getPath(new URI('http://www.nextflow.io/index.html'))

        when:
        def attrs = fsp.readHttpAttributes(path)
        then:
        attrs.lastModifiedTime() == null
        attrs.size() > 0
    }

    def "should read file attributes from FtpPath"() {
        given:
        def fsp = new FtpFileSystemProvider()
        def path = (XPath)fsp.getPath(new URI('ftp://ftp.ebi.ac.uk/robots.txt'))

        when:
        def attrs = fsp.readHttpAttributes(path)
        then:
        attrs.lastModifiedTime() == null
        attrs.size() == -1
    }
}
