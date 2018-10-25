/*
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

package nextflow.file.http

import spock.lang.Specification
import spock.lang.Unroll

/**
 * @author Emilio Palumbo
 * @author Paolo DiTommaso
 *
 */
class XPathTest extends Specification {

    def 'should validate equals and hashCode' () {

        when:
        def p1 = XPath.get('http://www.nextflow.io/a/b/c.txt')
        def p2 = XPath.get('http://www.nextflow.io/a/b/c.txt')
        def p3 = XPath.get('http://www.nextflow.io/z.txt')
        def p4 = XPath.get('http://www.google.com/a/b/c.txt')

        then:
        p1 == p2
        p1 != p3
        p1 != p4
        p1.equals(p2)
        !p1.equals(p3)
        !p1.equals(p4)
        p1.hashCode() == p2.hashCode()
        p1.hashCode() != p3.hashCode()
        p1.hashCode() != p4.hashCode()
    }


    def 'should convert to a string' () {
        expect:
        XPath.get('http://www.nextflow.io/abc/d.txt').toString()== '/abc/d.txt'
        XPath.get('http://www.nextflow.io/abc/d.txt').toUri() == new URI('http://www.nextflow.io/abc/d.txt')
        XPath.get('http://www.nextflow.io/abc/d.txt').toUri().toString()== 'http://www.nextflow.io/abc/d.txt'
    }

    def "should return url root"() {

        expect:
        XPath.get(origin).getRoot() == XPath.get(root)
        XPath.get(origin).getRoot().toString() == '/'
        XPath.get(origin).getRoot().toUri() == new URI(uri)

        where:
        origin                              | root                      | path      | uri
        'http://www.google.com/abc.txt'     | 'http://www.google.com/'  | '/'       | 'http://www.google.com/'
        'http://www.google.com/'            | 'http://www.google.com/'  | '/'       | 'http://www.google.com/'
        'http://www.google.com'             | 'http://www.google.com/'  | '/'       | 'http://www.google.com/'
    }


    def 'should return file name from url' () {

        expect:
        XPath.get(origin)?.getFileName() == XPath.get(fileName)
        XPath.get(origin)?.getFileName()?.toString() == fileName

        where:
        origin                          | fileName
        'http://nextflow.io/a/b/c.txt'  | 'c.txt'
        'http://nextflow.io/alpha'      | 'alpha'
        'http://nextflow.io/'           | null
        'http://nextflow.io'            | null
    }


    def 'should return if it is an absolute path'() {

        expect:
        XPath.get(origin).isAbsolute() == expected

        where:
        origin                          | expected
        'http://nextflow.io/a/b/c.txt'  | true
        'http://nextflow.io/a/b/'       | true
        'name.txt'                      | false

    }

    def 'should return parent path' () {

        expect:
        XPath.get(origin).getParent() == XPath.get(parent)

        where:
        origin                          | parent
        'http://nextflow.io/a/b/c.txt'  | 'http://nextflow.io/a/b/'
        'http://nextflow.io/a/'         | 'http://nextflow.io/'
        'http://nextflow.io/a'          | 'http://nextflow.io/'
        'http://nextflow.io/'           | null
    }


    @Unroll
    def 'should return name count #origin' () {

        expect:
        XPath.get(origin).getNameCount() == count

        where:
        origin                          | count
        'http://nextflow.io/a/b/c.txt'  | 3
        'http://nextflow.io/a/b/'       | 2
        'http://nextflow.io/a/b'        | 2
        'http://nextflow.io/a/'         | 1
        'http://nextflow.io/'           | 0
        'http://nextflow.io'            | 0
        'hello/world'                   | 2
        'hello'                         | 1

    }

    def 'should return name part by index' () {
        expect:
        XPath.get('http://nextflow.io/a/b/c.txt').getName(0) == XPath.get('a')
        XPath.get('http://nextflow.io/a/b/c.txt').getName(1) == XPath.get('b')
        XPath.get('http://nextflow.io/a/b/c.txt').getName(2) == XPath.get('c.txt')

        when:
        XPath.get('http://nextflow.io/a/b/c.txt').getName(3)
        then:
        thrown(IllegalArgumentException)

    }

    def 'should return a subpath' () {

        expect:
        XPath.get('http://nextflow.io/a/b/c/d.txt').subpath(0,3) == XPath.get('a/b/c' )
        XPath.get('http://nextflow.io/a/b/c/d.txt').subpath(1,3) == XPath.get('b/c' )
        XPath.get('http://nextflow.io/a/b/c/d.txt').subpath(3,4) == XPath.get('d.txt' )
    }

    @Unroll
    def 'should resolve a path: #base' () {

        expect:
        XPath.get(base).resolve(ext) == XPath.get(expected)
        XPath.get(base).resolve(XPath.get(ext)) == XPath.get(expected)

        where:
        base                        | ext                   | expected
        'http://nextflow.io/abc'    | 'd.txt'               | 'http://nextflow.io/abc/d.txt'
        'http://nextflow.io/abc/'   | 'd.txt'               | 'http://nextflow.io/abc/d.txt'
        'http://nextflow.io'        | 'file.txt'            | 'http://nextflow.io/file.txt'
        'http://nextflow.io/'       | 'file.txt'            | 'http://nextflow.io/file.txt'
        'alpha'                     | 'beta.txt'            | 'alpha/beta.txt'
        '/alpha'                    | 'beta.txt'            | '/alpha/beta.txt'
        'http://nextflow.io/abc/'   | '/z.txt'              | 'http://nextflow.io/z.txt'
        'http://nextflow.io/abc/'   | 'http://x.com/z.txt'  | 'http://x.com/z.txt'

    }

    @Unroll
    def 'should resolve sibling: #base' () {

        expect:
        XPath.get(base).resolveSibling(ext) == XPath.get(expected)
        XPath.get(base).resolveSibling(XPath.get(ext)) == XPath.get(expected)

        where:
        base                         | ext           | expected
        'http://nextflow.io/ab/c'    | 'd.txt'       | 'http://nextflow.io/ab/d.txt'
        'http://nextflow.io/ab/c/'   | 'd.txt'       | 'http://nextflow.io/ab/d.txt'
        'http://nextflow.io/ab'      | 'd.txt'       | 'http://nextflow.io/d.txt'
        'http://nextflow.io/ab/'     | 'd.txt'       | 'http://nextflow.io/d.txt'
        'http://nextflow.io/'        | 'd.txt'       | 'd.txt'
        'http://nextflow.io/a/b/c/'  | '/p/q.txt'    | 'http://nextflow.io/p/q.txt'
        'http://nextflow.io/abc/'    | 'http://x.com/z.txt'  | 'http://x.com/z.txt'
    }

    def 'should normalise a path: #base' () {
        expect:
        XPath.get(path).normalize() == XPath.get(expected)

        where:
        path                                | expected
        'http://nextflow.io/ab/c.txt'       | 'http://nextflow.io/ab/c.txt'
        'http://nextflow.io/ab/c/../d.txt'  | 'http://nextflow.io/ab/d.txt'
        'ab/c/../d.txt'                     | 'ab/d.txt'

    }

    def 'should iterate over a path' () {

        when:
        def itr = XPath.get('http://nextflow.io/a/b/c.txt').iterator()
        then:
        itr.hasNext()
        itr.next() == XPath.get('a')
        itr.hasNext()
        itr.next() == XPath.get('b')
        itr.hasNext()
        itr.next() == XPath.get('c.txt')
        !itr.hasNext()
        itr.next() == null
    }

}

