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
        appender.add('eng', 'Hello')
        appender.add('ita', 'Ciao')
        appender.add('eng', ' world!')
        appender.add('ita', '\nmondo\n!')
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
        appender.add('eng', 'Hello')
        appender.add('ita', 'Ciao')
        appender.add('eng', 'world!')
        appender.add('ita', 'mondo')
        appender.add('ita', '!')
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
        appender1.add('ENG', 'Hello')
        appender1.add('ITA', 'Ciao')
        appender1.add('ENG', 'world!')
        appender1.add('ITA', 'mondo')
        appender1.add('ITA', '!')

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
        appender2.add('ENG', 'Hello')
        appender2.add('ITA', 'Ciao')
        appender2.add('ENG', 'world!')
        appender2.add('ITA', 'mondo')
        appender2.add('ITA', '!')

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


    def testAppenderWithSeedClosure() {

        when:
        def file1 = Files.createTempFile('noname',null)
        file1.text = 'Inizio'

        def appender1 = new SimpleFileCollector()
        appender1.newLine = true
        appender1.seed = { key -> key == 'ENG' ? 'English' : 'Italian' }
        appender1.add('ENG', 'Hello')
        appender1.add('ITA', 'Ciao')
        appender1.add('ENG', 'world!')
        appender1.add('ITA', 'mondo')
        appender1.add('ITA', '!')

        then:
        appender1.get('ENG').text == 'English\nHello\nworld!\n'
        appender1.get('ITA').text == 'Italian\nCiao\nmondo\n!\n'
        appender1.size() == 2

        cleanup:
        file1?.delete()
        appender1?.close()

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
        appender.add('x', file1).add('x', '\n').add('x', file2)
        appender.add('y', file2).add('y', '\n').add('y', file3)


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
        appender.add('eng', 'Hello')
        appender.add('ita', 'Ciao')
        appender.add('eng', ' world!')
        appender.add('ita', '\nmondo\n!')
        appender.add('xxx', file1)
        then:
        appender.size() == 3
        appender.get('eng').text == 'Hello world!'
        appender.get('ita').text == 'Ciao\nmondo\n!'
        appender.get('xxx').text == 'file-content'
        file1.exists()

        when:
        def target = Files.createTempDirectory('new-dir')
        def list = appender.saveTo(target)
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

    def 'should append csv chunks with header' () {

        given:
        def appender = new SimpleFileCollector(keepHeader: true, skipLines: 1)

        when:
        appender.add('foo', 'COL1,COL2,COL3\naaa,bbb,ccc\nppp,qqq,rrr\n')
        appender.add('foo', 'COL1,COL2,COL3\nvvv,www,sss\nxxx,yyy,zzz\n')
        appender.add('foo', 'COL1,COL2,COL3\n111,222,333\n444,555,666\n')

        then:
        appender.size()==1
        appender.get('foo').text == '''
            COL1,COL2,COL3
            aaa,bbb,ccc
            ppp,qqq,rrr
            vvv,www,sss
            xxx,yyy,zzz
            111,222,333
            444,555,666
            '''
            .stripIndent().leftTrim()

    }

    def 'should append csv chunks with more headers' () {

        given:
        def appender = new SimpleFileCollector(keepHeader: true, skipLines: 3)

        when:
        appender.add('foo', '#\n#\nCOL1,COL2,COL3\naaa,bbb,ccc\nppp,qqq,rrr\n')
        appender.add('bar', '#\n#\nNUM1,NUM2,NUM3\n111,222,333\n444,555,666\n')
        appender.add('foo', '#\n#\nCOL1,COL2,COL3\nvvv,www,sss\nxxx,yyy,zzz\n')
        appender.add('bar', '#\n#\nNUM1,NUM2,NUM3\n777,888,999\n000,111,222\n')

        then:
        appender.size()==2

        appender.get('foo').text == '''
            #
            #
            COL1,COL2,COL3
            aaa,bbb,ccc
            ppp,qqq,rrr
            vvv,www,sss
            xxx,yyy,zzz
            '''
                .stripIndent().leftTrim()

        appender.get('bar').text == '''
            #
            #
            NUM1,NUM2,NUM3
            111,222,333
            444,555,666
            777,888,999
            000,111,222
            '''
                .stripIndent().leftTrim()

    }


}