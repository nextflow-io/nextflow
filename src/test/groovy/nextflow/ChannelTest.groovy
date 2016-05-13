/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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
import static java.nio.file.StandardWatchEventKinds.ENTRY_CREATE
import static java.nio.file.StandardWatchEventKinds.ENTRY_DELETE
import static java.nio.file.StandardWatchEventKinds.ENTRY_MODIFY

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import org.junit.Rule
import spock.lang.Specification
import test.TemporaryPath

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ChannelTest extends Specification {

    def setupSpec() {
        new Session()
    }


    def testSingleFile() {

        when:
        def channel = Channel.fromPath('/some/file.txt')
        then:
        channel.val == Paths.get('/some/file.txt')

        when:
        channel = Channel.fromPath('/some/f{i}le.txt')
        then:
        channel.val == Paths.get('/some/f{i}le.txt')

    }

    @Rule
    TemporaryPath tempDir = new TemporaryPath()

    def testGlobAlternative() {

        setup:
        def folder = tempDir.root
        def file1 = Files.createFile(folder.resolve('alpha.txt'))
        def file2 = Files.createFile(folder.resolve('beta.txt'))
        def file3 = Files.createFile(folder.resolve('gamma.txt'))
        def file4 = Files.createFile(folder.resolve('file4.txt'))
        def file5 = Files.createFile(folder.resolve('file5.txt'))
        def file6 = Files.createFile(folder.resolve('file66.txt'))

        when:
        def result = Channel
                        .fromPath("$folder/{alpha,gamma}.txt")
                        .toSortedList().getVal().collect { it -> it.name }
        then:
        result == [ 'alpha.txt', 'gamma.txt' ]

        when:
        result = Channel
                        .fromPath("$folder/file?.txt")
                        .toSortedList().getVal().collect { it -> it.name }
        then:
        result == [ 'file4.txt', 'file5.txt' ]

        when:
        result = Channel
                    .fromPath("$folder/file*.txt")
                    .toSortedList().getVal().collect { it -> it.name }
        then:
        result == [ 'file4.txt', 'file5.txt', 'file66.txt' ]

        when:
        result = Channel
                .fromPath("$folder/{alpha,}.txt")
                .toSortedList().getVal().collect { it -> it.name }
        then:
        result == [ 'alpha.txt']

    }


    def testGlobHiddenFiles() {

        setup:
        final folder = tempDir.root
        def file1 = Files.createFile(folder.resolve('.alpha.txt'))
        def file2 = Files.createFile(folder.resolve('.beta.txt'))
        def file3 = Files.createFile(folder.resolve('delta.txt'))
        def file4 = Files.createFile(folder.resolve('gamma.txt'))

        when:
        def result = Channel.fromPath("$folder/*").toSortedList().getVal()
        then:
        result == [file3, file4]

        when:
        result = Channel.fromPath("$folder/.*").toSortedList().getVal()
        then:
        result == [file1, file2]

        when:
        result = Channel.fromPath("$folder/{.*,*}", hidden: true).toSortedList().getVal()
        then:
        result == [file1, file2, file3, file4]

    }

    def testGlobFiles() {

        setup:
        def folder = tempDir.root
        def file1 = Files.createFile(folder.resolve('file1.txt'))
        def file2 = Files.createFile(folder.resolve('file2.txt'))
        def file3 = Files.createFile(folder.resolve('file3.txt'))
        def file4 = Files.createFile(folder.resolve('file4.log'))
        def sub1 = Files.createDirectories(folder.resolve('sub1'))
        def file5 = Files.createFile(sub1.resolve('file5.log'))
        def file6 = Files.createFile(sub1.resolve('file6.txt'))

        when:
        def result = Channel.fromPath("$folder/*.txt").toSortedList().getVal()
        then:
        result == [file1, file2, file3]

        when:
        def result2 = Channel.fromPath("$folder/**.txt").toSortedList().getVal()
        then:
        result2 == [file1, file2, file3, file6]

        when:
        def result3 = Channel.fromPath("$folder/sub1/**.log").toSortedList().getVal()
        then:
        result3 == [file5]

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

    def testFromPath() {

        setup:
        def folder = tempDir.root
        def file1 = Files.createFile(folder.resolve('file1.txt'))
        def file2 = Files.createFile(folder.resolve('file2.txt'))
        def file3 = Files.createFile(folder.resolve('file3.log'))
        def sub1 = Files.createDirectories(folder.resolve('sub1'))
        def file5 = Files.createFile(sub1.resolve('file5.log'))

        when:
        List<Path> result = Channel
                                .fromPath( folder.toAbsolutePath().toString() + '/*.txt' )
                                .toSortedList().getVal().collect { it.name }
        then:
        result == [ 'file1.txt', 'file2.txt' ]

        when:
        result = Channel
                    .fromPath( folder.toAbsolutePath().toString() + '/*' )
                    .toSortedList().getVal().collect { it.name }
        then:
        result == [ 'file1.txt', 'file2.txt', 'file3.log' ]

        when:
        result = Channel
                    .fromPath( folder.toAbsolutePath().toString() + '/*', type: 'file' )
                    .toSortedList().getVal().collect { it.name }
        then:
        result == [ 'file1.txt', 'file2.txt', 'file3.log' ]

        when:
        result = Channel
                    .fromPath( folder.toAbsolutePath().toString() + '/*', type: 'dir' )
                    .toSortedList().getVal().collect { it.name }
        then:
        result == ['sub1']

        when:
        result = Channel
                    .fromPath( folder.toAbsolutePath().toString() + '/*', type: 'any' )
                    .toSortedList().getVal().collect { it.name }
        then:
        result == [ 'file1.txt', 'file2.txt', 'file3.log', 'sub1' ]

        when:
        result = Channel
                    .fromPath( folder.toAbsolutePath().toString() + '/**', type: 'file' )
                    .toSortedList() .getVal() .collect { it.name }
        then:
        result == [ 'file1.txt', 'file2.txt', 'file3.log', 'file5.log' ]

        when:
        def result2 = Channel
                    .fromPath( folder.toAbsolutePath().toString() + '/**', type: 'file', maxDepth: 0 )
                    .toSortedList() .getVal() .collect { it.name }
        then:
        result2 == ['file1.txt', 'file2.txt', 'file3.log' ]

        when:
        def result3 = Channel
                        .fromPath (folder.toAbsolutePath().toString() + '/{file1.txt,sub1/file5.log}')
                        .toSortedList() .val .collect { it.name }
        then:
        result3 == ['file1.txt','file5.log']

        when:
        Channel.fromPath( folder.toAbsolutePath().toString() + '/*', xx: 'any' )
        then:
        thrown( IllegalArgumentException )

        when:
        Channel.fromPath( folder.toAbsolutePath().toString() + '/*', type: 'ciao' )
        then:
        thrown( IllegalArgumentException )

    }


    def testFromPathWithLinks() {

        setup:
        def folder = tempDir.root
        def file1 = Files.createFile(folder.resolve('file1.txt'))
        def file2 = Files.createFile(folder.resolve('file2.txt'))
        def sub1 = Files.createDirectories(folder.resolve('sub_1'))
        def file3 = Files.createFile(sub1.resolve('file3.txt'))
        def file4 = Files.createFile(sub1.resolve('file4.txt'))
        Files.createSymbolicLink(folder.resolve('link_to_sub1'), sub1 )

        // -- by default traverse symlinks
        when:
        def result = Channel.fromPath( folder.toAbsolutePath().toString() + '/**/*.txt' ).toSortedList({it.name}).getVal().collect { it.getName() }
        then:
        result == ['file3.txt','file3.txt','file4.txt','file4.txt']

        // -- switch off symlinks traversing
        when:
        def result2 = Channel.fromPath( folder.toAbsolutePath().toString() + '/**/*.txt', followLinks: false ).toSortedList({it.name}).getVal().collect { it.getName() }
        then:
        result2 == ['file3.txt','file4.txt']

    }


    def testFromPathHidden() {

        setup:
        def folder = tempDir.root
        Files.createFile(folder.resolve('file1.txt'))
        Files.createFile(folder.resolve('file2.txt'))
        Files.createFile(folder.resolve('.file_hidden.txt'))

        // -- by default no hidden
        when:
        def result = Channel.fromPath( folder.toAbsolutePath().toString() + '/*.txt' ).toSortedList({it.name}).getVal().collect { it.getName() }
        then:
        result == ['file1.txt','file2.txt']

        when:
        result = Channel.fromPath( folder.toAbsolutePath().toString() + '/.*.txt' ).toSortedList({it.name}).getVal().collect { it.getName() }
        then:
        result == ['.file_hidden.txt']

        when:
        result = Channel.fromPath( folder.toAbsolutePath().toString() + '/*.txt', hidden: true ).toSortedList({it.name}).getVal().collect { it.getName() }
        then:
        result == ['.file_hidden.txt', 'file1.txt','file2.txt']


    }


}