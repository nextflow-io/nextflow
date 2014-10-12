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

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
import java.nio.file.Files

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SimpleFileCollectorTest extends Specification {


    def testAppenderWithString() {

        when:
        def appender = new SimpleFileCollector()
        appender.append('eng', 'Hello')
        appender.append('ita', 'Ciao')
        appender.append('eng', ' world!')
        appender.append('ita', '\nmondo\n!')
        then:
        appender.size() == 2
        appender.get('eng').text == 'Hello world!'
        appender.get('ita').text == 'Ciao\nmondo\n!'
        appender.getFiles().size() == 2
        appender.getFiles() *. name .sort() == ['eng','ita']

        cleanup:
        appender?.close()

    }

    def testAppenderNewLine() {

        when:
        def appender = new SimpleFileCollector()
        appender.newLine = true
        appender.append('eng', 'Hello')
        appender.append('ita', 'Ciao')
        appender.append('eng', 'world!')
        appender.append('ita', 'mondo')
        appender.append('ita', '!')
        then:
        appender.size() == 2
        appender.get('eng').text == 'Hello\nworld!\n'
        appender.get('ita').text == 'Ciao\nmondo\n!\n'
        appender.getFiles().size() == 2
        appender.getFiles() *. name .sort() == ['eng','ita']

        cleanup:
        appender?.close()

    }

    def testAppenderWithSeed() {

        when:
        def file1 = Files.createTempFile('noname',null)
        file1.text = 'Inizio'

        def appender1 = new SimpleFileCollector()
        appender1.newLine = true
        appender1.seed = [ENG: 'Begin', ITA: file1]
        appender1.append('ENG', 'Hello')
        appender1.append('ITA', 'Ciao')
        appender1.append('ENG', 'world!')
        appender1.append('ITA', 'mondo')
        appender1.append('ITA', '!')

        then:
        appender1.get('ENG').text == 'Begin\nHello\nworld!\n'
        appender1.get('ITA').text == 'Inizio\nCiao\nmondo\n!\n'
        appender1.size() == 2

        when:
        def file2 = Files.createTempFile('noname',null)
        file2.text = 'same file'

        def appender2 = new SimpleFileCollector()
        appender2.newLine = true
        appender2.seed = file2
        appender2.append('ENG', 'Hello')
        appender2.append('ITA', 'Ciao')
        appender2.append('ENG', 'world!')
        appender2.append('ITA', 'mondo')
        appender2.append('ITA', '!')

        then:
        appender2.get('ENG').text == 'same file\nHello\nworld!\n'
        appender2.get('ITA').text == 'same file\nCiao\nmondo\n!\n'
        appender2.size() == 2

        cleanup:
        file1?.delete()
        file2?.delete()
        appender1?.close()
        appender2?.close()

    }


    def testAppenderWithFile() {

        given:
        def file1 = Files.createTempFile('file1',null)
        file1.text = 'alpha\nbeta'

        def file2 = Files.createTempFile('file2',null)
        file2.text = 'Hello\nworld'

        def file3 = Files.createTempFile('file3',null)
        file3.text = 'xyz'

        when:
        def appender = new SimpleFileCollector()
        appender.append('x', file1).append('x', '\n').append('x', file2)
        appender.append('y', file2).append('y', '\n').append('y', file3)


        then:
        appender.size() == 2
        appender.get('x').text == 'alpha\nbeta\nHello\nworld'
        appender.get('y').text == 'Hello\nworld\nxyz'

        cleanup:
        appender?.close()
        file1?.delete()
        file2?.delete()
        file3?.delete()

    }

    def testMove() {
        when:
        def file1 = Files.createTempFile('testFile',null)
        file1.text = 'file-content'

        def appender = new SimpleFileCollector()
        appender.append('eng', 'Hello')
        appender.append('ita', 'Ciao')
        appender.append('eng', ' world!')
        appender.append('ita', '\nmondo\n!')
        appender.append('xxx', file1)
        then:
        appender.size() == 3
        appender.get('eng').text == 'Hello world!'
        appender.get('ita').text == 'Ciao\nmondo\n!'
        appender.get('xxx').text == 'file-content'
        file1.exists()

        when:
        def target = Files.createTempDirectory('new-dir')
        def list = appender.moveFiles(target)
        then:
        list.size() == 3
        list *. name .sort() == ['eng','ita','xxx']
        appender.size() == 0
        file1.exists()

        when:
        appender.close()
        then:
        file1.exists()

        cleanup:
        appender?.close()
        target?.deleteDir()
        file1?.delete()
    }



}