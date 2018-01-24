/*
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
