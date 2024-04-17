/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.extension

import java.nio.file.Path
import java.nio.file.Paths
import java.text.SimpleDateFormat
import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneId

import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Ignore
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
        def defLocale = Locale.getDefault(Locale.Category.FORMAT)
        def useLocale = new Locale.Builder().setLanguage(locale).build()
        Locale.setDefault(Locale.Category.FORMAT, useLocale)

        when:
        def formatter = new SimpleDateFormat(fmt, Locale.ENGLISH)
        formatter.setTimeZone(TimeZone.getTimeZone(tz))
        then:
        Bolts.format(now, fmt, tz) == expected
        formatter.format(now) == expected

        cleanup:
        Locale.setDefault(Locale.Category.FORMAT, defLocale)

        where:
        tz      | fmt                       | locale  | expected
        'UTC'   | 'dd-MM-yyyy HH:mm'        | 'en'    | '14-12-2017 21:12'
        'CET'   | 'dd-MM-yyyy HH:mm'        | 'en'    | '14-12-2017 22:12'
        'UTC'   | 'dd-MMM-yyyy HH:mm:ss'    | 'en'    | '14-Dec-2017 21:12:27'
        'CST'   | 'dd-MM-yyyy HH:mm'        | 'en'    | '14-12-2017 15:12'
        'CST'   | 'dd-MMM-yyyy HH:mm:ss'    | 'en'    | '14-Dec-2017 15:12:27'
        'UTC'   | 'dd-MMM-yyyy HH:mm:ss'    | 'es'    | '14-Dec-2017 21:12:27'
        'CST'   | 'dd-MMM-yyyy HH:mm:ss'    | 'es'    | '14-Dec-2017 15:12:27'
    }

    @Ignore("we dont need to test java functionalities")
    def 'should format offset datetime' () {
        given:
        def now = OffsetDateTime.ofInstant(Instant.ofEpochMilli(1513285947928), ZoneId.of('CET'))
        expect:
        now.format(FMT) == EXPECTED

        where:
        FMT                     | EXPECTED
        'dd-MM-yyyy HH:mm'      | '14-12-2017 22:12'
        'dd-MMM-yyyy HH:mm:ss'  | '14-Dec-2017 22:12:27'

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

    def testAsURL() {
        expect: 
        'http://foo.com' as URL == new URL('http://foo.com')
        'http://foo.com/some/file.txt' as URL == new URL('http://foo.com/some/file.txt')
    }

    def testAsURI() {
        expect:
        'http://foo.com' as URI == URI.create('http://foo.com')
        'http://foo.com/some/file.txt' as URI == URI.create('http://foo.com/some/file.txt')
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

    def 'should convert num to memory' () {
        expect:
        2_048.toMemory() instanceof MemoryUnit
        1_048_576.toMemory() == MemoryUnit.of('1 MB')
    }

    def 'should convert num to duration' () {
        expect:
        1_000.toDuration() instanceof Duration
        60_000.toDuration() == Duration.of('1 min')
    }

    def 'should get or create a value' () {

        given:
        def map = [X: 1]

        when:
        def val = map.getOrCreate('X', 9)
        then:
        val == 1
        map.size() == 1

        when:
        val = map.getOrCreate('Y', 9)
        then:
        val == 9
        map.size() == 2
        map.X == 1
        map.Y == 9

    }

    def 'should get or create a value with closure' () {

        given:
        def map = [X: 1]

        when:
        def val = map.getOrCreate('X') { return 9 }
        then:
        val == 1
        map.size() == 1

        when:
        val = map.getOrCreate('Y') { return 9 }
        then:
        val == 9
        map.size() == 2
        map.X == 1
        map.Y == 9

    }

    def 'should redact secrets ' () {

        expect:
        Bolts.redact('foo') == '...'
        Bolts.redact('12345') == '...'
        Bolts.redact('123456') == '1...'
        Bolts.redact('123456789') == '1234...'
        Bolts.redact('12345678901234567890') == '12345...'
        and:
        Bolts.redact('foo', 3) == '...'
        Bolts.redact('12345', 3) == '12...'
        and:
        Bolts.redact('foo', 3, 'xx') == 'xx'
        Bolts.redact('12345', 3, 'xx') == '12xx'
    }

    def 'should deep clone obj' () {
        given:
        Map map = [foo: 1, bar: [x: 2, y: 3]]

        when:
        def copy = Bolts.deepClone0((Serializable)map)
        then:
        copy == map

        when:
        copy.bar.x = 20
        then:
        copy.bar.x == 20
        map.bar.x == 2
    }

    def 'should deep clone map' () {
        given:
        Map map = [foo: 1, bar: [x: 2, y: 3]]

        when:
        def copy = Bolts.deepClone((Map)map)
        then:
        copy == map

        when:
        copy.bar.x = 20
        then:
        copy.bar.x == 20
        map.bar.x == 2
    }

    def 'should deep merge map and overwrite values'() {
        given:
        def origMap = [foo: 1, bar: [x: 2, y: 3]]
        and:
        def newMap = [bar: [x: 4, z: 5]]

        when:
        def merge = Bolts.deepMerge(origMap, newMap)

        then:
        merge.foo == 1 
        merge.bar.x == 4
        merge.bar.y == 3
        merge.bar.z == 5
    }

    def 'should deep merge config object and overwrite values'() {
        given:
        def bar = new ConfigObject(); bar.putAll([x: 2, y: 3])
        def cfg = new ConfigObject(); cfg.putAll( [foo: 1, bar: bar] )
        and:
        def map = [bar: [x: 4, z: 5]]

        when:
        def result = Bolts.deepMerge(cfg, map)

        then:
        result.foo == 1
        result.bar.x == 4
        result.bar.y == 3
        result.bar.z == 5
        and:
        result instanceof ConfigObject
        result.bar instanceof ConfigObject
        and:
        !result.is( cfg )
    }

}
