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

package nextflow

import org.junit.Rule
import spock.lang.Specification
import spock.lang.Timeout

import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.TimeUnit

import groovyx.gpars.dataflow.DataflowQueue
import test.TemporaryPath
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(value = 1, unit = TimeUnit.MINUTES)
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

    def testFromPathWithListOfPatterns() {

        setup:
        def folder = tempDir.root.toAbsolutePath()
        folder.resolve('file1.txt').text = 'Hello'
        folder.resolve('file2.txt').text = 'Hola'
        folder.resolve('file3.fq').text = 'Ciao'

        when:
        def result = Channel.fromPath( ["$folder/*.txt", "$folder/*.fq"] ).toSortedList({it.name}).getVal().collect { it.getName() }
        then:
        result == [ 'file1.txt', 'file2.txt', 'file3.fq' ]

        when:
        result = Channel.fromPath( ["$folder/file1.txt", "$folder/file2.txt", "$folder/file3.fq"] ).toSortedList({it.name}).getVal().collect { it.getName() }
        then:
        result == [ 'file1.txt', 'file2.txt', 'file3.fq' ]

        when:
        result = Channel.fromPath( ["$folder/*"] ).toSortedList({it.name}).getVal().collect { it.getName() }
        then:
        result == [ 'file1.txt', 'file2.txt', 'file3.fq' ]
    }

    def 'should check file exists' () {

        given:
        def folder = tempDir.root.toAbsolutePath()

        def file1 = folder.resolve('file1.txt')
        def file2 = folder.resolve('file2.txt')
        def file3 = folder.resolve('file3.txt')

        // only file1 and file3 exist
        file1.text = 'foo'
        file3.text = 'bar'

        when:
        def ch = Channel.fromPath(file1.toString())
        then:
        ch.getVal() == file1
        ch.getVal() == Channel.STOP

        when:
        ch = Channel.fromPath([file1.toString(), file2.toString()])
        then:
        ch.getVal() == file1
        ch.getVal() == file2
        ch.getVal() == Channel.STOP

        when:
        ch = Channel.fromPath(file1.toString(), checkIfExists: true)
        then:
        ch.getVal() == file1
        ch.getVal() == Channel.STOP

        when:
        def session = new Session()
        Channel.fromPath(file2.toString(), checkIfExists: true)
        Channel.fromPath0Future.get()
        then:
        session == Global.session
        session.aborted
        session.error instanceof NoSuchFileException
        session.error.message == file2.toString()

        when:
        session = new Session()
        Channel.fromPath([file1, file2, file3], checkIfExists: true)
        Channel.fromPath0Future.get()
        then:
        session.aborted
        session.error instanceof NoSuchFileException
        session.error.message == file2.toString()

        when:
        session = new Session()
        Channel.fromPath('http://google.com/foo.txt', checkIfExists: true)
        Channel.fromPath0Future.get()
        then:
        session.aborted
        session.error instanceof NoSuchFileException
        session.error.message == 'http://google.com/foo.txt'
    }

    def 'should check if pattern exists' () {

        given:
        def folder = Files.createTempDirectory('test'); folder.deleteDir()
        def session = new Session()
        when:
        Channel.fromPath("$folder/*.txt", checkIfExists: true)
        Channel.fromPath0Future.get()
        then:
        session.isAborted()
        session.getError() instanceof IllegalArgumentException

    }

    def 'should check if pattern exists for multiple patterns' () {

        given:
        def session = new Session()
        def folder = tempDir.root.toAbsolutePath()
        def file1 = folder.resolve('file1.txt'); file1.text = 'foo'

        when:
        def result = Channel.fromPath("$folder/*.txt", checkIfExists: true)
        then:
        result.getVal() instanceof Path
        !session.terminated

        when:
        session = new Session()
        Channel.fromPath(["$folder/*.txt", 'foo/*.fq'], checkIfExists: true)
        Channel.fromPath0Future.get()
        then:
        session.isAborted()
        session.getError() instanceof IllegalArgumentException

        cleanup:
        folder?.deleteDir()
    }

    def 'should check if pattern exists for filePairs' () {

        when:
        def session = new Session()
        Channel.fromFilePairs('foo/*.txt', checkIfExists: true)
        Channel.fromPath0Future.get()
        then:
        session.isAborted()
        session.getError() instanceof IllegalArgumentException

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
        Channel.readPrefix(Paths.get('/some/path/abc_1.fa'), 'abc_1.fa') == 'abc_1'

        when:
        Channel.readPrefix(Paths.get('/some/path/abc_1.fa'), 'foo')
        then:
        thrown(IllegalArgumentException)
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

    def 'should group files with the same prefix using a custom grouping' () {

        setup:
        def folder = TestHelper.createInMemTempDir()
        def a1 = Files.createFile(folder.resolve('hello_1.fa'))
        def a2 = Files.createFile(folder.resolve('hello_2.fa'))
        def b1 = Files.createFile(folder.resolve('hola_1.fa'))
        def b2 = Files.createFile(folder.resolve('hola_2.fa'))
        def c1 = Files.createFile(folder.resolve('ciao_1.fa'))
        def c2 = Files.createFile(folder.resolve('ciao_2.fa'))

        when:
        def grouping = { Path file -> file.name.substring(0,1) }
        def pairs = Channel.fromFilePairs(folder.resolve("*_{1,2}.*"), grouping, size:-1)
        then:
        pairs.val == ['c', [c1, c2]]
        pairs.val == ['h', [a1, a2, b1, b2]]
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

    def 'should group file with a collection of patterns' () {

        given:
        def folder = TestHelper.createInMemTempDir()
        def a1 = Files.createFile(folder.resolve('alpha_1.fa'))
        def a2 = Files.createFile(folder.resolve('alpha_2.fa'))
        def b1 = Files.createFile(folder.resolve('beta_1.fa'))
        def b2 = Files.createFile(folder.resolve('beta_2.fa'))

        def d1 = Files.createFile(folder.resolve('delta_1.fq'))
        def d2 = Files.createFile(folder.resolve('delta_2.fq'))
        def g1 = Files.createFile(folder.resolve('gamma_1.fq'))
        def g2 = Files.createFile(folder.resolve('gamma_2.fq'))

        when:
        def pairs = Channel.fromFilePairs([folder.resolve("*_{1,2}.fa"), folder.resolve("$folder/*_{1,2}.fq")])
        then:
        pairs.val == ['alpha', [a1, a2]]
        pairs.val == ['beta', [b1, b2]]
        pairs.val == ['delta', [d1, d2]]
        pairs.val == ['gamma', [g1, g2]]
        pairs.val == Channel.STOP

    }

    def 'should return singleton group' () {
        given:
        def folder = TestHelper.createInMemTempDir()
        def a1 = Files.createFile(folder.resolve('alpha_1.fa'))
        def a2 = Files.createFile(folder.resolve('alpha_2.fa'))

        when:
        def files = Channel.fromFilePairs(folder.resolve('*.fa'), size:1)
        then:
        files.val == ['alpha_1', [a1]]
        files.val == ['alpha_2', [a2]]
        files.val == Channel.STOP
    }

    def 'should return singleton' () {
        given:
        def folder = TestHelper.createInMemTempDir()
        def a1 = Files.createFile(folder.resolve('alpha_1.fa'))

        when:
        def files = Channel.fromFilePairs(a1, size:1)
        then:
        files.val == ['alpha_1', [a1]]
        files.val == Channel.STOP
    }

    def 'should use size one by default' () {
        given:
        def folder = TestHelper.createInMemTempDir()
        def a1 = Files.createFile(folder.resolve('alpha_1.fa'))

        when:
        def files = Channel.fromFilePairs(a1)
        then:
        files.val == ['alpha_1', [a1]]
        files.val == Channel.STOP
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
        given:
        def folder = Files.createTempDirectory('test')
        when:
        def result = Channel.watchPath("$folder/foo/*")
        then:
        result.val == Channel.STOP
        cleanup:
        folder?.deleteDir()
    }

}
