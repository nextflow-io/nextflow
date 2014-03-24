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

package nextflow
import java.nio.file.Files
import java.nio.file.Paths
import static java.nio.file.StandardWatchEventKinds.*

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ChannelTest extends Specification {

    def testGetPathAndPattern () {

        expect:
        Channel.getFolderAndPattern( '/some/file/name.txt' ) == ['/some/file/', 'name.txt']
        Channel.getFolderAndPattern( '/some/file/*.txt' ) == ['/some/file/', '*.txt']
        Channel.getFolderAndPattern( '/some/file/*' ) == ['/some/file/', '*']
        Channel.getFolderAndPattern( '/some/file/' ) == ['/some/file/', '']
        Channel.getFolderAndPattern( 'path/filename.txt' ) == ['path/', 'filename.txt']
        Channel.getFolderAndPattern( 'filename.txt' ) == ['./', 'filename.txt']
        Channel.getFolderAndPattern( './file.txt' ) == ['./', 'file.txt']

        Channel.getFolderAndPattern( '/some/file/**/*.txt' ) == ['/some/file/', '**/*.txt']

        Channel.getFolderAndPattern( 'dxfs:///some/file/**/*.txt' ) == ['dxfs:///some/file/', '**/*.txt']
        Channel.getFolderAndPattern( 'dxfs://some/file/**/*.txt' ) == ['dxfs://some/file/', '**/*.txt']
        Channel.getFolderAndPattern( 'dxfs://*.txt' ) == ['dxfs://./', '*.txt']
        Channel.getFolderAndPattern( 'dxfs:///*.txt' ) == ['dxfs:///', '*.txt']
        Channel.getFolderAndPattern( 'dxfs:///**/*.txt' ) == ['dxfs:///', '**/*.txt']
    }


    def testSingleFile() {

        when:
        def channel = Channel.path('/some/file.txt')
        then:
        channel.val == Paths.get('/some/file.txt')

        when:
        channel = Channel.path('/some/f{i}le.txt')
        then:
        channel.val == Paths.get('/some/f{i}le.txt')

    }


    def testGlobAlternative() {

        setup:
        def folder = Files.createTempDirectory('testFiles')
        def file1 = Files.createFile(folder.resolve('alpha.txt'))
        def file2 = Files.createFile(folder.resolve('beta.txt'))
        def file3 = Files.createFile(folder.resolve('gamma.txt'))
        def file4 = Files.createFile(folder.resolve('file4.txt'))
        def file5 = Files.createFile(folder.resolve('file5.txt'))
        def file6 = Files.createFile(folder.resolve('file66.txt'))

        when:
        def channel = Channel.path("$folder/{alpha,gamma}.txt")
        then:
        channel.val == folder.resolve('alpha.txt')
        channel.val == folder.resolve('gamma.txt')
        channel.val == Channel.STOP

        when:
        channel = Channel.path("$folder/file?.txt")
        then:
        channel.val == folder.resolve('file4.txt')
        channel.val == folder.resolve('file5.txt')
        channel.val == Channel.STOP

        when:
        channel = Channel.path("$folder/file*.txt")
        then:
        channel.val == folder.resolve('file4.txt')
        channel.val == folder.resolve('file5.txt')
        channel.val == folder.resolve('file66.txt')
        channel.val == Channel.STOP

        cleanup:
        folder.deleteDir()

    }


    def testGlobHiddenFiles() {

        setup:
        def folder = Files.createTempDirectory('testFiles')
        def file1 = Files.createFile(folder.resolve('.alpha.txt'))
        def file2 = Files.createFile(folder.resolve('.beta.txt'))
        def file3 = Files.createFile(folder.resolve('delta.txt'))
        def file4 = Files.createFile(folder.resolve('gamma.txt'))

        when:
        def channel = Channel.path("$folder/*")
        then:
        channel.val == folder.resolve('delta.txt')
        channel.val == folder.resolve('gamma.txt')
        channel.val == Channel.STOP

        when:
        channel = Channel.path("$folder/.*")
        then:
        channel.val == folder.resolve('.alpha.txt')
        channel.val == folder.resolve('.beta.txt')
        channel.val == Channel.STOP


        when:
        channel = Channel.path("$folder/{.*,*}")
        then:
        channel.val == folder.resolve('.alpha.txt')
        channel.val == folder.resolve('.beta.txt')
        channel.val == folder.resolve('delta.txt')
        channel.val == folder.resolve('gamma.txt')
        channel.val == Channel.STOP

        cleanup:
        folder.deleteDir()

    }

    def testGlobFiles() {

        setup:
        def folder = Files.createTempDirectory('testFiles')
        def file1 = Files.createFile(folder.resolve('file1.txt'))
        def file2 = Files.createFile(folder.resolve('file2.txt'))
        def file3 = Files.createFile(folder.resolve('file3.txt'))
        def file4 = Files.createFile(folder.resolve('file4.log'))
        def sub1 = Files.createDirectories(folder.resolve('sub1'))
        def file5 = Files.createFile(sub1.resolve('file5.log'))
        def file6 = Files.createFile(sub1.resolve('file6.txt'))

        when:
        def channel = Channel.path("$folder/*.txt")
        then:
        channel.val == file1
        channel.val == file2
        channel.val == file3
        channel.val == Channel.STOP


        when:
        def channel2 = Channel.path("$folder/**.txt")
        then:
        channel2.val.toString() == file1.toString()
        channel2.val.toString() == file2.toString()
        channel2.val.toString() == file3.toString()
        channel2.val.toString() == file6.toString()
        channel2.val == Channel.STOP

        when:
        def channel3 = Channel.path("$folder/sub1/**.log")
        then:
        channel3.val.toString() == file5.toString()
        channel3.val == Channel.STOP

        cleanup:
        folder.deleteDir()

    }

    def testChopLines() {

        when:
        def result = Channel.readLines('a\nbb\nccc') { it.reverse() }
        then:
        result.val == 'a'
        result.val == 'bb'
        result.val == 'ccc'
        result.val == Channel.STOP

    }

    def testChopFasta() {

        setup:
        String fasta =
            /
            >1aboA
            NLFVALYDFVASGD
            NTLSITKGEKL
            >1ycsB
            RVLGYNHNGEWCEA
            QTKNGQGWVPS
            /.stripIndent()


        when:
        def result = Channel.readFasta(fasta, record:[id:true])
        then:
        result.val .id == '1aboA'
        result.val .id == '1ycsB'
        result.val == Channel.STOP

    }

    def testStringEvents() {

        when:
        Channel.stringToWatchEvents('xxx')
        then:
        thrown(IllegalArgumentException)

        expect:
        Channel.stringToWatchEvents() == [ ENTRY_CREATE ]
        Channel.stringToWatchEvents('create,delete') == [ENTRY_CREATE, ENTRY_DELETE]
        Channel.stringToWatchEvents('Create , MODIFY ') == [ENTRY_CREATE, ENTRY_MODIFY]

    }


}