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

package nextflow.extension
import java.nio.file.Path
import java.nio.file.Paths
import java.text.SimpleDateFormat

import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BoltsTest extends Specification {

    @Unroll
    def 'should format a date' () {

        given:
        def now = new Date(1513285947928)

        when:
        def formatter = new SimpleDateFormat(fmt)
        formatter.setTimeZone(TimeZone.getTimeZone(tz))
        then:
        Bolts.format(now, fmt, tz) == expected
        formatter.format(now) == expected

        where:
        tz      | fmt                       | expected
        'UTC'   | 'dd-MM-yyyy HH:mm'        | '14-12-2017 21:12'
        'CET'   | 'dd-MM-yyyy HH:mm'        | '14-12-2017 22:12'
        'UTC'   | 'dd-MMM-yyyy HH:mm:ss'    | '14-Dec-2017 21:12:27'
        'CST'   | 'dd-MM-yyyy HH:mm'        | '14-12-2017 15:12'
        'CST'   | 'dd-MMM-yyyy HH:mm:ss'    | '14-Dec-2017 15:12:27'

    }

    def testTrimDotZero() {

        expect:
        ''.trimDotZero() == ''
        '.0'.trimDotZero() == ''
        '1.0'.trimDotZero() == '1'
        '123.000'.trimDotZero() == '123'

    }


    def 'test leftTrim' () {

        expect:
        '  hola hello  '.leftTrim() == 'hola hello  '
        '\n\n hola hello\n'.leftTrim() == 'hola hello\n'

    }

    def 'test rightTrim' () {

        expect:
        '  hola hello  '.rightTrim() == '  hola hello'
        '\n\nhola hello\n\n'.rightTrim() == '\n\nhola hello'

    }


    def testAsDuration() {

        setup:
        def x = 3;

        expect:
        2_000 as Duration == Duration.of('2 second')
        '1s' as Duration == Duration.of('1 second')
        "$x min" as Duration == Duration.of('3 min')

    }

    def testAsPath() {

        setup:
        def x = 'Hello'

        expect:
        null as Path == null
        'file.txt' as Path == Paths.get('file.txt')
        '/some/path/file.txt' as Path == Paths.get('/some/path/file.txt')
        "name.fa" as Path == Paths.get('name.fa')
        "/some/path/${x}.txt" as Path == Paths.get('/some/path/Hello.txt')

        new File('/path/to/file') as Path == Paths.get('/path/to/file')


    }

    def testAsMemoryUnit() {
        given:
        def x = 5

        expect:
        1024 as MemoryUnit == MemoryUnit.of('1 KB')
        '10 GB' as MemoryUnit == MemoryUnit.of('10 GB')
        "$x MB" as MemoryUnit == MemoryUnit.of('5 MB')
    }

    def testConfigToMap  () {

        setup:
        def config = new ConfigSlurper().parse( 'task {field1=1; field2="two"}; env { x = 99 }; list { q = [1,2,3]  }' )

        when:
        def map = Bolts.toMap(config)
        map.env.PATH = '/local/bin'

        then:
        !(map instanceof ConfigObject)
        map.task.field1 == 1
        map.task.field2 == "two"
        map.env.x == 99
        map.env.y == null
        map.env.PATH == '/local/bin'

        map.list.q == [1,2,3]

    }


    def testBestMatches() {

        expect:
        ['hello','hola','ciao'].bestMatches('ciao') == ['ciao']
        ['hello','hola','ciao'].bestMatches('halo') == ['hello','hola']
        ['hello','hola','ciao'].bestMatches('none') == []

    }


    def testNavigate() {

        when:
        def map = [ a: 1, b:2 ]
        then:
        map.navigate('a') == 1
        map.navigate('b') == 2
        map.navigate('b') { it == 2 ? 4 : 0 }  == 4
        map.navigate('c') == null
        map.navigate('x.y.z') == null


        when:
        map = [ a: 1, b:2, c: [x:1, y: [delta: 'd', gamma: 'g']] ]
        then:
        map.navigate('a') == 1
        map.navigate('b') == 2
        map.navigate('c.x') == 1
        map.navigate('c.y.delta') == 'd'
        map.navigate('c.y.gamma') == 'g'
        map.navigate('c.y.zeta') == null
        map.navigate('c.p.q') == null
        map.navigate('c.y.gamma') { it == 'g' ? 'omega' : null }== 'omega'

        when:
        map = [trace: [:]]
        then:
        map.navigate('trace.x') == null
        map.navigate('trace.x') { return 1 } == null


    }

    def testIsCamelCase() {

        expect:
        !'hello'.isCamelCase()
        'helloWorld'.isCamelCase()
        !'hello_world'.isCamelCase()
        !'HELLOworld'.isCamelCase()
        !'Helloworld'.isCamelCase()
        'HelloWorld'.isCamelCase()

    }

    def 'should return the nearest match' () {

        given:
        def salut = ['hello','hola','halo','hey','ciao','bonjour']

        expect:
        salut.closest('hola') == ['hola']
        salut.closest('hol') == ['hola']
        salut.closest('cioa') == ['ciao']
        salut.closest('helo') == ['hello', 'halo']
    }

    def 'should check lower case strings' () {
        expect:
        'abc'.isLowerCase()
        !'AAA'.isLowerCase()
        '1a11a1a'.isLowerCase()
    }

    def 'should check upper case strings' () {
        expect:
        'ABCD'.isUpperCase()
        !'aaaa'.isUpperCase()
        '1A11A1A'.isUpperCase()
    }

    def 'should check all lower case strings' () {
        expect:
        'abc'.isAllLowerCase()
        !'AAA'.isAllLowerCase()
        !'1a11a1a'.isAllLowerCase()
    }

    def 'should check all upper case strings' () {
        expect:
        'ABCD'.isAllUpperCase()
        !'aaaa'.isAllUpperCase()
        !'1A11A1A'.isAllUpperCase()
    }

    def 'should convert a map to a list of pairs' () {
        expect:
        [a:1, b:2, c:3].pairs() == [['a',1], ['b',2], ['c',3]]
        [a:1, b:[2,4], c:[3,9]].pairs() == [['a',1], ['b',[2,4]], ['c',[3,9]]]
        [a:1, b:[2,4], c:[3,9]].pairs(flat:true) == [['a',1], ['b',2], ['b',4], ['c',3], ['c',9]]
    }

}
