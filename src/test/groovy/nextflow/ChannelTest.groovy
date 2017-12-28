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

package nextflow

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import groovyx.gpars.dataflow.DataflowQueue
import org.junit.Rule
import spock.lang.Specification
import test.TemporaryPath
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ChannelTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def testFrom() {
        given:
        DataflowQueue result

        when:
        result = Channel.from('hola')
        then:
        result.val == 'hola'
        result.val == Channel.STOP

        when:
        result = Channel.from('alpha','delta')
        then:
        result.val == 'alpha'
        result.val == 'delta'
        result.val == Channel.STOP

        when:
        result = Channel.from(['alpha','delta'])
        then:
        result.val == 'alpha'
        result.val == 'delta'
        result.val == Channel.STOP

        when:
        result = Channel.from([])
        then:
        result.val == Channel.STOP
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
        result == ['alpha.txt']

        when:
        result = Channel
                .fromPath("$folder/{,beta}.txt")
                .toSortedList().getVal().collect { it -> it.name }
        then:
        result == ['beta.txt']

        when:
        result = Channel
                .fromPath("$folder/alpha.txt{,}")
                .toSortedList().getVal().collect { it -> it.name }
        then:
        result == ['alpha.txt']
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

    def testEscapeGlob() {

        setup:
        def folder = tempDir.root
        def file1 = Files.createFile(folder.resolve('file1.txt'))
        def file2 = Files.createFile(folder.resolve('file*.txt'))
        def file3 = Files.createFile(folder.resolve('file?.txt'))
        def sub1 = Files.createDirectories(folder.resolve('sub[a-b]'))
        def file5 = Files.createFile(sub1.resolve('file5.log'))
        def file6 = Files.createFile(sub1.resolve('file6.txt'))

        when:
        def result = Channel.fromPath("$folder/file\\*.txt").toSortedList().getVal()
        then:
        result == [file2]

        when:
        result = Channel.fromPath("$folder/file*.txt", glob: false).toSortedList().getVal()
        then:
        result == [file2]

        when:
        result = Channel.fromPath("$folder/sub\\[a-b\\]/file*").toSortedList().getVal()
        then:
        result == [file5,file6]

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


    def testFromPathS3() {

        when:
        Channel.fromPath('s3://bucket/some/data.txt')
        then:
        noExceptionThrown()
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


    def 'should return files prefix' () {

        expect:
        Channel.readPrefix(Paths.get('/some/path/abc1.fa'), '*1.fa') == 'abc'
        Channel.readPrefix(Paths.get('/some/path/abc_1.fa'), '*1.fa') == 'abc'
        Channel.readPrefix(Paths.get('/some/path/abc_1.fa'), '*_1.fa') == 'abc'
        Channel.readPrefix(Paths.get('/some/path/abc_1.fa'), '*_{1,2}.fa') == 'abc'
        Channel.readPrefix(Paths.get('/some/path/abc_1.fa'), '*_[1-2].fa') == 'abc'
        Channel.readPrefix(Paths.get('/some/path/abc_1.fa'), '**_{1,2}.*') == 'abc'
        Channel.readPrefix(Paths.get('/some/path/abc_1_trimmed.fa'), '*_[1-2]_*.fa') == 'abc'
        Channel.readPrefix(Paths.get('/some/path/abc_1_trimmed.fa'), '*??_[1-2]_*.fa') == 'abc'
        Channel.readPrefix(Paths.get('/some/path/abc_1.fa'), 'abc_{1,2}.fa') == 'abc'
        Channel.readPrefix(Paths.get('/some/path/abc_1.fa'), 'abc_[1-9].fa') == 'abc'
        Channel.readPrefix(Paths.get('/some/path/foo_abc_1.fa'), 'foo_*_{1,2}.fa') == 'foo_abc'
    }

    def 'should group files with the same prefix' () {

        setup:
        def folder = TestHelper.createInMemTempDir()
        def a1 = Files.createFile(folder.resolve('alpha_1.fa'))
        def a2 = Files.createFile(folder.resolve('alpha_2.fa'))
        def b1 = Files.createFile(folder.resolve('beta_1.fa'))
        def b2 = Files.createFile(folder.resolve('beta_2.fa'))
        def d1 = Files.createFile(folder.resolve('delta_1.fa'))
        def d2 = Files.createFile(folder.resolve('delta_2.fa'))

        when:
        def pairs = Channel.fromFilePairs(folder.resolve("*_{1,2}.*"))
        then:
        pairs.val == ['alpha', [a1, a2]]
        pairs.val == ['beta', [b1, b2]]
        pairs.val == ['delta', [d1, d2]]
        pairs.val == Channel.STOP

        when:
        pairs = Channel.fromFilePairs(folder.resolve("*_{1,2}.fa") , flat: true)
        then:
        pairs.val == ['alpha', a1, a2]
        pairs.val == ['beta', b1, b2]
        pairs.val == ['delta', d1, d2]
        pairs.val == Channel.STOP
    }

    def 'should group files with the same prefix and setting size' () {

        setup:
        def folder = TestHelper.createInMemTempDir()
        def a1 = Files.createFile(folder.resolve('alpha_1.fa'))
        def a2 = Files.createFile(folder.resolve('alpha_2.fa'))
        def a3 = Files.createFile(folder.resolve('alpha_3.fa'))

        def b1 = Files.createFile(folder.resolve('beta_1.fa'))
        def b2 = Files.createFile(folder.resolve('beta_2.fa'))
        def b3 = Files.createFile(folder.resolve('beta_3.fa'))

        def d1 = Files.createFile(folder.resolve('delta_1.fa'))
        def d2 = Files.createFile(folder.resolve('delta_2.fa'))
        def d3 = Files.createFile(folder.resolve('delta_3.fa'))
        def d4 = Files.createFile(folder.resolve('delta_4.fa'))

        when:
        // default size == 2
        def pairs = Channel.fromFilePairs(folder.resolve("*_{1,2,3}.fa"))
        then:
        pairs.val == ['alpha', [a1, a2]]
        pairs.val == ['beta', [b1, b2]]
        pairs.val == ['delta', [d1, d2]]
        pairs.val == Channel.STOP

        when:
        pairs = Channel.fromFilePairs(folder.resolve("*_{1,2,3,4}.fa"), size: 3)
        then:
        pairs.val == ['alpha', [a1, a2, a3]]
        pairs.val == ['beta', [b1, b2, b3]]
        pairs.val == ['delta', [d1, d2, d3]]
        pairs.val == Channel.STOP

        when:
        pairs = Channel.fromFilePairs(folder.resolve("*_{1,2,3,4}.fa"), size: -1)
        then:
        pairs.val == ['alpha', [a1, a2, a3]]
        pairs.val == ['beta', [b1, b2, b3]]
        pairs.val == ['delta', [d1, d2, d3, d4]]
        pairs.val == Channel.STOP

    }

    def 'should watch and emit a file' () {
        given:
        def folder = Files.createTempDirectory('test')

        when:
        def result = Channel.watchPath("$folder/")
        sleep 500
        Files.createFile(folder.resolve('hello.txt'))
        then:
        result.val == folder.resolve('hello.txt')

        when:
        result = Channel.watchPath(folder.toString())
        sleep 500
        Files.createFile(folder.resolve('ciao.txt'))
        then:
        result.val == folder.resolve('ciao.txt')

        cleanup:
        folder?.deleteDir()
    }

    def 'should return an empty channel when watching a missing path' () {

        when:
        def result = Channel.watchPath("foo/*")
        then:
        result.val == Channel.STOP
    }

}
