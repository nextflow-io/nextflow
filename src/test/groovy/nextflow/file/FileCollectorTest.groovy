/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.file
import java.nio.file.Files

import org.mapdb.DataInput2
import org.mapdb.DataOutput2
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FileCollectorTest extends Specification {


    def testSerializer () {

        given:
        final HELLO_WORLD = 'Hello world!'
        def buffer = new ByteArrayOutputStream()
        def dataOut = new DataOutputStream(buffer)
        new FileCollector.ObjectSerializer().serialize(dataOut, HELLO_WORLD)
        dataOut.close()

        when:
        def bytes = buffer.toByteArray()
        def dataIn = new DataInputStream(new ByteArrayInputStream(bytes))
        then:
        new FileCollector.ObjectSerializer().deserialize(dataIn, HELLO_WORLD.size()) == HELLO_WORLD

    }


    def testSerializer2() {

        given:
        def dataOut = new DataOutput2(new byte[1024])
        new FileCollector.ObjectSerializer().serialize(dataOut, 'Hello world!')
        dataOut.close()

        when:
        def buffer = dataOut.copyBytes()
        def dataIn = new DataInput2(buffer)
        then:
        new FileCollector.ObjectSerializer().deserialize(dataIn, buffer.length) == 'Hello world!'

    }

    def testAppenderWithString0() {

        given:
        def target = Files.createTempDirectory('test')

        def appender = new FileCollector()
        appender.append('eng', 'Hello')

        when:
        def files = appender.moveFiles(target)

        then:
        files.size() == 1
        target.resolve('eng').text == 'Hello'

        cleanup:
        appender?.close()
        target?.deleteDir()

    }


    def testAppenderWithString() {

        given:
        def target = Files.createTempDirectory('test')

        def appender = new FileCollector()
        appender.append('eng', 'Hello')
        appender.append('ita', 'Ciao')
        appender.append('eng', ' world!')
        appender.append('ita', '\nmondo\n!')

        when:
        def files = appender.moveFiles(target)

        then:
        files.size() == 2
        files *. name .sort() == ['eng','ita']
        target.resolve('eng').text == 'Hello world!'
        target.resolve('ita').text == 'Ciao\nmondo\n!'

        cleanup:
        appender?.close()
        target?.deleteDir()

    }

    def testAppenderNewLine() {

        given:
        def target = Files.createTempDirectory('test')
        def appender = new FileCollector()
        appender.newLine = true
        appender.append('eng', 'Hello')
        appender.append('ita', 'Ciao')
        appender.append('eng', 'world!')
        appender.append('ita', 'mondo')
        appender.append('ita', '!')

        when:
        def files = appender.moveFiles(target)

        then:
        files *. name .sort() == ['eng','ita']
        target.resolve('eng').text == 'Hello\nworld!\n'
        target.resolve('ita').text == 'Ciao\nmondo\n!\n'

        cleanup:
        appender?.close()
        target?.deleteDir()

    }

    def testAppenderWithSeed() {

        given:
        def target = Files.createTempDirectory('test')

        when:
        def file1 = Files.createTempFile('noname',null)
        file1.text = 'Inizio'

        def appender1 = new FileCollector()
        appender1.newLine = true
        appender1.seed = [ENG: 'Begin', ITA: file1]
        appender1.append('ENG', 'Hello')
        appender1.append('ITA', 'Ciao')
        appender1.append('ENG', 'world!')
        appender1.append('ITA', 'mondo')
        appender1.append('ITA', '!')

        // move files
        def files = appender1.moveFiles(target)

        then:
        files.size()==2
        target.resolve('ENG').text == 'Begin\nHello\nworld!\n'
        target.resolve('ITA').text == 'Inizio\nCiao\nmondo\n!\n'


        when:
        target.deleteDir()
        def file2 = Files.createTempFile('noname',null)
        file2.text = 'same file'

        def appender2 = new FileCollector()
        appender2.newLine = true
        appender2.seed = file2
        appender2.append('ENG', 'Hello')
        appender2.append('ITA', 'Ciao')
        appender2.append('ENG', 'world!')
        appender2.append('ITA', 'mondo')
        appender2.append('ITA', '!')

        files = appender2.moveFiles(target)

        then:
        files.size() == 2
        target.resolve('ENG').text == 'same file\nHello\nworld!\n'
        target.resolve('ITA').text == 'same file\nCiao\nmondo\n!\n'

        cleanup:
        file1?.delete()
        file2?.delete()
        target?.deleteDir()

        appender1?.close()
        appender2?.close()

    }


    def testAppenderWithFile() {

        given:
        def target = Files.createTempDirectory('test')

        def file1 = Files.createTempFile('file1',null)
        file1.text = 'alpha\nbeta'

        def file2 = Files.createTempFile('file2',null)
        file2.text = 'Hello\nworld'

        def file3 = Files.createTempFile('file3',null)
        file3.text = 'xyz'

        when:
        def appender = new FileCollector()
        appender.append('x', file1).append('x', '\n').append('x', file2)
        appender.append('y', file2).append('y', '\n').append('y', file3)
        def files = appender.moveFiles(target)

        then:
        files.size() == 2
        target.resolve('x').text == 'alpha\nbeta\nHello\nworld'
        target.resolve('y').text == 'Hello\nworld\nxyz'

        cleanup:
        appender?.close()
        file1?.delete()
        file2?.delete()
        file3?.delete()
        target?.deleteDir()

    }


    def testAppenderWithSort() {

        given:
        def target = Files.createTempDirectory('test')

        def appender = new FileCollector()
        appender.sort = { String it -> it.length() }
        appender.newLine = true

        appender.append('file2', 'bbbbbbb')
        appender.append('file1', 'cccccc')
        appender.append('file2', 'aaaaa')
        appender.append('file1', 'zzzz')
        appender.append('file2', 'pp')
        appender.append('file1', 'q')

        when:
        def files = appender.moveFiles(target)

        then:
        files.size() == 2
        files *. name .sort() == ['file1','file2']
        target.resolve('file1').text == 'q\nzzzz\ncccccc\n'
        target.resolve('file2').text == 'pp\naaaaa\nbbbbbbb\n'

        cleanup:
        appender?.close()
        target?.deleteDir()

    }

}
