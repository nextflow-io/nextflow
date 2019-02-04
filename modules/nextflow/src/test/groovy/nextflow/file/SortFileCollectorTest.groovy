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

package nextflow.file
import java.nio.file.Files

import nextflow.exception.AbortOperationException
import org.iq80.leveldb.DB
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SortFileCollectorTest extends Specification {

    def 'test append strings and save files'() {

        given:
        def folder = Files.createTempDirectory('test')

        when:
        def collector = new SortFileCollector()
        collector.add('alpha', 'BBB')
        collector.add('delta', '222')
        collector.add('alpha', 'AAA')
        collector.add('delta', '111')
        collector.saveTo(folder)
        /*
         * No sorting has been specified
         * Entries are appended in the order they have been added
         */
        then:
        folder.list().size() == 2
        folder.resolve('alpha').text == 'BBBAAA'
        folder.resolve('delta').text == '222111'

        cleanup:
        collector?.close()
        folder?.deleteDir()

    }

    def 'test append strings and save sorted files'() {

        given:
        def folder = Files.createTempDirectory('test')

        when:
        def collector = new SortFileCollector()
        collector.sort = { it -> it }
        collector.add('alpha', 'BBB')
        collector.add('delta', '222')
        collector.add('alpha', 'AAA')
        collector.add('delta', '111')
        collector.saveTo(folder)
        then:
        folder.list().size() == 2
        folder.resolve('alpha').text == 'AAABBB'
        folder.resolve('delta').text == '111222'

        cleanup:
        collector?.close()
        folder?.deleteDir()

    }


    def 'test append strings with new-line separator and seed values'() {

        given:
        def folder = Files.createTempDirectory('test')

        when:
        def collector = new SortFileCollector()
        collector.sort = { it -> it }
        collector.newLine = true
        collector.seed = [alpha: '000', delta: '111', gamma: '222']
        collector.add('alpha', 'BBB')
        collector.add('delta', 'qqq')
        collector.add('alpha', 'ZZZ')
        collector.add('delta', 'ttt')
        collector.add('gamma', 'yyy')
        collector.add('gamma', 'zzz')
        collector.add('gamma', 'xxx')
        collector.add('delta', 'ppp')
        collector.add('alpha', 'AAA')
        collector.saveTo(folder)
        then:
        folder.list().size() == 3
        folder.resolve('alpha').text == '000\nAAA\nBBB\nZZZ\n'
        folder.resolve('delta').text == '111\nppp\nqqq\nttt\n'
        folder.resolve('gamma').text == '222\nxxx\nyyy\nzzz\n'

        cleanup:
        collector?.close()
        folder?.deleteDir()

    }


    def 'test append strings with new-line separator and seed closure'() {

        given:
        def folder = Files.createTempDirectory('test')

        when:
        def collector = new SortFileCollector()
        collector.sort = { it -> it }
        collector.newLine = true
        collector.seed = { key -> key=='alpha' ? '000' : (key=='delta' ? '111' : '222') }
        collector.add('alpha', 'BBB')
        collector.add('delta', 'qqq')
        collector.add('alpha', 'ZZZ')
        collector.add('delta', 'ttt')
        collector.add('gamma', 'yyy')
        collector.add('gamma', 'zzz')
        collector.add('gamma', 'xxx')
        collector.add('delta', 'ppp')
        collector.add('alpha', 'AAA')
        collector.saveTo(folder)
        then:
        folder.list().size() == 3
        folder.resolve('alpha').text == '000\nAAA\nBBB\nZZZ\n'
        folder.resolve('delta').text == '111\nppp\nqqq\nttt\n'
        folder.resolve('gamma').text == '222\nxxx\nyyy\nzzz\n'

        cleanup:
        collector?.close()
        folder?.deleteDir()

    }

    def 'test sort file collect properties'() {

        given:
        def folder = Files.createTempDirectory('test')
        def identity = { it -> it }

        when:
        def collector = new SortFileCollector()
        collector.sort = identity
        collector.newLine = true
        collector.seed = [alpha: '000', delta: '111', gamma: '222']
        collector.sliceMaxItems = 100
        collector.sliceMaxSize = 20_000
        collector.deleteTempFilesOnClose = false

        collector.add('alpha', 'BBB')
        collector.add('delta', 'qqq')
        collector.add('alpha', 'ZZZ')
        collector.add('delta', 'ttt')
        collector.add('gamma', 'yyy')
        collector.add('gamma', 'zzz')
        collector.add('gamma', 'xxx')
        collector.add('delta', 'ppp')
        collector.add('alpha', 'AAA')
        collector.saveTo(folder)
        then:
        folder.list().size() == 3
        folder.resolve('alpha').text == '000\nAAA\nBBB\nZZZ\n'
        folder.resolve('delta').text == '111\nppp\nqqq\nttt\n'
        folder.resolve('gamma').text == '222\nxxx\nyyy\nzzz\n'
        collector.index.sliceMaxItems == 100
        collector.index.sliceMaxSize == 20_000
        collector.index.getTempDir() == collector.getTempDir().resolve("index")
        collector.index.deleteTempFilesOnClose == collector.deleteTempFilesOnClose
        (collector.index.comparator as SortFileCollector.IndexSort).sort == identity

        cleanup:
        collector?.close()
        folder?.deleteDir()
        collector?.getTempDir()?.deleteDir()

    }


    def 'test sort file collector with closure'() {

        given:
        def folder = Files.createTempDirectory('test')

        when:
        def collector = new SortFileCollector()
        collector.sort = { it.size() }
        collector.newLine = true
        collector.add('x', 'AAAA')
        collector.add('x', 'AAAAAA')
        collector.add('x', 'A')
        collector.add('x', 'AA')
        collector.saveTo(folder)
        /*
         * No sorting has been specified
         * Entries are appended in the order they have been added
         */
        then:
        folder.list().size() == 1
        folder.resolve('x').text == 'A\nAA\nAAAA\nAAAAAA\n'

        cleanup:
        collector?.close()
        folder?.deleteDir()

    }

    def 'test sort file collector with comparator'() {

        given:
        def folder = Files.createTempDirectory('test')

        when:
        def collector = new SortFileCollector()
        collector.sort = { o1, o2 -> o2 <=> o1 } as Comparator
        collector.newLine = true
        collector.add('x', 'A')
        collector.add('x', 'B')
        collector.add('x', 'C')
        collector.add('x', 'D')
        collector.saveTo(folder)
        /*
         * No sorting has been specified
         * Entries are appended in the order they have been added
         */
        then:
        folder.list().size() == 1
        folder.resolve('x').text == 'D\nC\nB\nA\n'

        cleanup:
        collector?.close()
        folder?.deleteDir()

    }

    def 'test create sort comparator' () {

        def collector
        SortFileCollector.IndexSort criteria

        when:
        collector = new SortFileCollector()
        collector.store = Mock(DB)
        criteria = collector.createSortComparator()
        then:
        criteria.sort == null
        criteria.comparator == null

        when:
        def closure = { -> }
        collector = new SortFileCollector()
        collector.store = Mock(DB)
        collector.sort = closure
        criteria = collector.createSortComparator()
        then:
        criteria.sort == closure
        criteria.comparator == null


        when:
        def comp = { a, b -> a <=> b } as Comparator
        collector = new SortFileCollector()
        collector.store = Mock(DB)
        collector.sort = comp
        criteria = collector.createSortComparator()
        then:
        criteria.sort == null
        criteria.comparator == comp

        when:
        collector = new SortFileCollector()
        collector.store = Mock(DB)
        collector.sort = 'any'
        then:
        thrown(AbortOperationException)

    }

    def 'should collect file and keep headers' () {

        given:
        def folder = Files.createTempDirectory('test')
        def appender = new SortFileCollector(keepHeader: true, skipLines: 1)

        when:
        appender.add('foo', 'COL1,COL2,COL3\naaa,bbb,ccc\nppp,qqq,rrr\n')
        appender.add('bar', 'NUM1,NUM2,NUM3\n111,222,333\n444,555,666\n')
        appender.add('foo', 'COL1,COL2,COL3\nvvv,www,sss\nxxx,yyy,zzz\n')
        appender.add('bar', 'NUM1,NUM2,NUM3\n777,888,999\n000,111,222\n')
        appender.saveTo(folder)

        then:
        folder.list().size()==2

        folder.resolve('foo').text == '''
            COL1,COL2,COL3
            aaa,bbb,ccc
            ppp,qqq,rrr
            vvv,www,sss
            xxx,yyy,zzz
            '''
                .stripIndent().leftTrim()

        folder.resolve('bar').text == '''
            NUM1,NUM2,NUM3
            111,222,333
            444,555,666
            777,888,999
            000,111,222
            '''
                .stripIndent().leftTrim()


        cleanup:
        appender?.close()
        folder?.deleteDir()

    }

}
