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

import java.nio.file.Files
import java.nio.file.Path

import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Channel
import nextflow.Session
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CollectFileOperatorTest extends Specification {

    @Shared
    def Session session

    def setupSpec() {
        session = new Session(cacheable: true)
        session.workDir = Files.createTempDirectory('nxf-test')
    }

    def cleanupSpec() {
        session.workDir.deleteDir()
    }

    def testCollectFileString() {

        when:
        def result = Channel
                .from('alpha','beta','gamma')
                .collectFile { it == 'beta' ? ['file2', it.reverse() ] : ['file1',it] }
                .toSortedList { it.name }

        List<Path> list = result.val

        then:
        list[0].name == 'file1'
        list[0].text == 'alphagamma'

        list[1].name == 'file2'
        list[1].text == 'ateb'

    }


    def testCollectFileWithFiles() {

        given:
        def file1 = Files.createTempDirectory('temp').resolve('A')
        file1.deleteOnExit()
        file1.text = 'alpha\nbeta'

        def file2 = Files.createTempDirectory('temp').resolve('B')
        file2.deleteOnExit()
        file2.text = 'Hello\nworld'

        def file3 = Files.createTempDirectory('temp').resolve('A')
        file3.deleteOnExit()
        file3.text = 'xyz'

        when:
        def list = Channel
                .from(file1,file2,file3)
                .collectFile(sort:'index')
                .toSortedList { it.name }
                .getVal() as List<Path>

        then:
        list[0].name == 'A'
        list[0].text == 'alpha\nbetaxyz'

        list[1].name == 'B'
        list[1].text == 'Hello\nworld'


        when:
        list = Channel
                .from(file1,file2,file3)
                .collectFile(sort:'index', newLine:true)
                .toSortedList { it.name }
                .getVal() as List<Path>

        then:
        list[0].name == 'A'
        list[0].text == 'alpha\nbeta\nxyz\n'

        list[1].name == 'B'
        list[1].text == 'Hello\nworld\n'


    }

    def testCollectManyFiles() {


        when:
        def list = Channel
                .from('Hola', 'Ciao', 'Hello', 'Bonjour', 'Halo')
                .collectFile(sort:'index') { item -> [ "${item[0]}.txt", item + '\n' ] }
                .toList()
                .getVal()
                .sort { it.name }

        then:
        list[0].name == 'B.txt'
        list[0].text == 'Bonjour\n'
        list[1].text == 'Ciao\n'
        list[1].name == 'C.txt'
        list[2].name == 'H.txt'
        list[2].text == 'Hola\nHello\nHalo\n'

    }


    def testCollectFileWithStrings() {

        when:
        def result = Channel
                .from('alpha', 'beta', 'gamma')
                .collectFile(name: 'hello.txt', newLine: true, sort:'index')

        def file = result.val

        then:
        result.val == Channel.STOP
        file.name == 'hello.txt'
        file.text == 'alpha\nbeta\ngamma\n'
    }

    def testCollectFileWithDefaultName() {

        when:
        def result = Channel
                .from('alpha', 'beta', 'gamma')
                .collectFile(newLine: true, sort:'index')

        def file = result.val

        then:
        result.val == Channel.STOP
        file.name.startsWith('collect')
        file.text == 'alpha\nbeta\ngamma\n'
    }

    def testCollectFileAndSortWithClosure() {

        when:
        def result = Channel
                .from('delta', 'beta', 'gamma','alpha')
                .collectFile(newLine: true, sort:{ it -> it })

        def file = result.val

        then:
        result.val == Channel.STOP
        file.name.startsWith('collect')
        file.text == 'alpha\nbeta\ndelta\ngamma\n'
    }

    def testCollectFileAndSortWithComparator() {

        when:
        def result = Channel
                .from('delta', 'beta', 'gamma','alpha')
                .collectFile(newLine: true, sort:{ a,b -> b<=>a } as Comparator)

        def file = result.val

        then:
        result.val == Channel.STOP
        file.name.startsWith('collect')
        file.text == 'gamma\ndelta\nbeta\nalpha\n'
    }


    def 'should collect file and skip header line' () {

        given:
        def file1 = Files.createTempDirectory('temp').resolve('A')
        file1.deleteOnExit()
        file1.text = 'HEADER\nalpha\nbeta\n'

        def file2 = Files.createTempDirectory('temp').resolve('B')
        file2.deleteOnExit()
        file2.text = 'HEADER\nHello\nworld\n'

        def file3 = Files.createTempDirectory('temp').resolve('A')
        file3.deleteOnExit()
        file3.text = 'HEADER\nxxx\nyyy\nzzz\n'


        when:
        def files = Channel
                .from(file1,file2,file3)
                .collectFile(skip:1, sort: 'index')
                .toList()
                .getVal()

        def result = [:]; files.each{ result[it.name]=it }
        then:
        result.A.name == 'A'
        result.A.text == 'alpha\nbeta\nxxx\nyyy\nzzz\n'

        result.B.name == 'B'
        result.B.text == 'Hello\nworld\n'


        when:
        files = Channel
                .from(file1,file2,file3)
                .collectFile(skip:2, sort: 'index')
                .toList()
                .getVal()

        result = [:]; files.each{ result[it.name]=it }
        then:
        result.A.name == 'A'
        result.A.text == 'beta\nyyy\nzzz\n'

        result.B.name == 'B'
        result.B.text == 'world\n'


        when:
        files = Channel
                .from(file1,file2,file3)
                .collectFile(skip:3, sort: 'index')
                .toList()
                .getVal()

        result = [:]; files.each{ result[it.name]=it }
        then:
        result.A.name == 'A'
        result.A.text == 'zzz\n'

        result.B.name == 'B'
        result.B.text == ''


        when:
        files = Channel
                .from(file1,file2,file3)
                .collectFile(skip:10, sort: 'index')
                .toList()
                .getVal()

        result = [:]; files.each{ result[it.name]=it }
        then:
        result.A.name == 'A'
        result.A.text == ''

        result.B.name == 'B'
        result.B.text == ''

    }

    def 'should collect file and keep header line' () {

        given:
        def file1 = Files.createTempDirectory('temp').resolve('A')
        file1.deleteOnExit()
        file1.text = 'HEADER\nalpha\nbeta\n'

        def file2 = Files.createTempDirectory('temp').resolve('B')
        file2.deleteOnExit()
        file2.text = '## HEAD ##\nHello\nworld\n'

        def file3 = Files.createTempDirectory('temp').resolve('A')
        file3.deleteOnExit()
        file3.text = 'HEADER\nxxx\nyyy\nzzz\n'


        when:
        def files = Channel
                .from(file1,file2,file3)
                .collectFile(keepHeader:true, sort: 'index')
                .toList()
                .getVal()

        def result = [:]; files.each{ result[it.name]=it }
        then:
        result.A.name == 'A'
        result.A.text == 'HEADER\nalpha\nbeta\nxxx\nyyy\nzzz\n'

        result.B.name == 'B'
        result.B.text == '## HEAD ##\nHello\nworld\n'

    }

    def 'check invalid options' () {
        given:
        CollectFileOp op

        when:
        op = new CollectFileOp(Mock(DataflowReadChannel), [seed: 'foo'])
        then:
        op.collector.seed == 'foo'

        when:
        op = new CollectFileOp(Mock(DataflowReadChannel), [skip: 20])
        then:
        op.collector.skipLines == 20

        when:
        op = new CollectFileOp(Mock(DataflowReadChannel), [keepHeader: true])
        then:
        op.collector.keepHeader
        op.collector.skipLines == 1

        when:
        op = new CollectFileOp(Mock(DataflowReadChannel), [keepHeader: true, skip: 10])
        then:
        op.collector.keepHeader
        op.collector.skipLines == 10

        when:
        new CollectFileOp(Mock(DataflowReadChannel), [seed: 'foo', keepHeader: true])
        then:
        thrown(IllegalArgumentException)

    }
}
