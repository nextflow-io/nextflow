/*
 * Copyright 2020-2022, Seqera Labs
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

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.ProcessStageException
import spock.lang.Ignore
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class FilePorterTest extends Specification {


    def 'should get the max retries value' () {

        when:
        def session = new Session()
        def porter = new FilePorter(session)
        then:
        porter.maxRetries == 3

        when:
        session = new Session([filePorter: [maxRetries: 88]])
        porter = new FilePorter(session)
        then:
        porter.maxRetries == 88

    }


    def 'should copy foreign files' () {

        given:
        def folder = Files.createTempDirectory('test')
        def session = new Session(workDir: folder)

        def foreign1 = TestHelper.createInMemTempFile('hola.txt', 'hola mundo!')
        def foreign2 = TestHelper.createInMemTempFile('ciao.txt', 'ciao mondo!')

        when:
        def porter = new FilePorter(session)
        def batch = porter.newBatch(folder)
        def file1 = batch.addToForeign(foreign1)
        def file2 = batch.addToForeign(foreign2)
        then:
        file1 != foreign1
        file1.name == foreign1.name
        and:
        file2 != foreign2
        file2.name == foreign2.name

        when:
        porter.transfer(batch)
        then:
        file1.name == 'hola.txt'
        file1.text == 'hola mundo!'
        file1.fileSystem == FileSystems.default
        and:
        file2.name == 'ciao.txt'
        file2.text == 'ciao mondo!'
        file2.fileSystem == FileSystems.default

        cleanup:
        folder?.deleteDir()
    }

    @Ignore
    def 'should not hang' () {
        given:
        def RANGE = (0..19)
        System.setProperty('filePorter.debugDelay', '5000')
        def folder = Files.createTempDirectory('test')
        def session = new Session(workDir: folder, threadPool: [FileTransfer: [maxThreads: 2]])

        List<Path> foreign1 = RANGE.collect { TestHelper.createInMemTempFile("hola-${it}.txt", "hola mundo ${it}") }
        List<Path> foreign2 = RANGE.collect { TestHelper.createInMemTempFile("ciao-${it}.txt", "ciao mondo ${it}") }
        List<Path> foreign3 = RANGE.collect { TestHelper.createInMemTempFile("helo-${it}.txt", "helo mondo ${it}") }

        and:
        def porter = new FilePorter(session)
        def batch1 = porter.newBatch(folder)
        foreign1.each { batch1.addToForeign(it) }

        and:
        def batch2 = porter.newBatch(folder)
        foreign2.each { batch2.addToForeign(it) }

        and:
        def batch3 = porter.newBatch(folder)
        foreign3.each { batch3.addToForeign(it) }

        when:
        def t1= Thread.start { porter.transfer(batch1) }
        and: sleep(2000)
        def t2= Thread.start { porter.transfer(batch2) }
        and: sleep(2000)
        def t3= Thread.start { porter.transfer(batch3) }

        then:
        noExceptionThrown()
        and:
        t1.join()
        t2.join()
        t3.join()

        cleanup:
        folder?.deleteDir()
    }



    static class ErrorStage extends FilePorter.FileTransfer {

        ErrorStage(Path path, Path stagePath, int maxRetries) {
            super(path, stagePath, maxRetries)
        }

        @Override
        void run() throws Exception {
            throw new ProcessStageException('Cannot stage gile')
        }
    }

    def 'should submit actions' () {

        given:
        def foreign1 = TestHelper.createInMemTempFile('hola.txt', 'hola mundo!')
        def foreign2 = TestHelper.createInMemTempFile('ciao.txt', 'ciao mondo!')

        def folder = Files.createTempDirectory('test')
        def session = new Session(workDir: folder)
        FilePorter porter = Spy(FilePorter, constructorArgs:[session])
        def files = [ foreign1, foreign2 ]

        when:
        def result = porter.submitStagingActions(files, folder)
        then:
        result.size() == 2
        result[0].result.done
        result[1].result.done

        when:
        porter.submitStagingActions([Paths.get('/missing/file')], folder)
        then:
        thrown(ProcessStageException)

        cleanup:
        folder?.deleteDir()
    }




    def 'should check batch size and truth' () {
        given:
        def STAGE = Files.createTempDirectory('test')
        def FTP_FILE = 'ftp://host.com/file.txt' as Path
        def MEM_FILE = TestHelper.createInMemTempFile('bar.txt')
        def session = Mock(Session)
        session.config >> [:]
        def porter = new FilePorter(session)

        when:
        def batch = porter.newBatch(STAGE)
        then:
        !batch
        batch.size() == 0

        when:
        def result1 = batch.addToForeign(MEM_FILE)
        then:
        result1.scheme == 'file'
        batch.size() == 1
        batch

        when:
        def result2 = batch.addToForeign(FTP_FILE)
        then:
        result2.scheme == 'file'
        batch.size() == 2
        batch

        cleanup:
        STAGE?.deleteDir()

    }

    def 'should return stage path' () {

        given:
        def STAGE = Files.createTempDirectory('test')
        def FTP_FILE1 = 'ftp://host.com/file1.txt' as Path
        def FTP_FILE2 = 'ftp://host.com/file2.txt' as Path

        when:
        def stage1 = FilePorter.getCachePathFor(FTP_FILE1, STAGE)
        then:
        stage1.toString().startsWith( STAGE.toString() )

        when:
        def stage2 = FilePorter.getCachePathFor(FTP_FILE2, STAGE)
        then:
        stage2.toString().startsWith( STAGE.toString() )
        stage1 != stage2

        when:
        STAGE.resolve('foo.txt').text = 'ciao' // <-- add a file to alter the content of the dir
        and:
        def newStage1 = FilePorter.getCachePathFor(FTP_FILE1, STAGE)
        then:
        stage1 == newStage1

        cleanup:
        STAGE?.deleteDir()
    }

    def 'should create transfer paths' () {
        given:
        def sess = Mock(Session) {
            getConfig() >> [:]
        }
        def folder = Files.createTempDirectory('test')
        def local = folder.resolve('hola.text'); local.text = 'Hola'
        and:
        def foreign = TestHelper.createInMemTempDir()

        and:
        def porter = new FilePorter(sess)
        
        when:
        def transfer1 = porter.createFileTransfer(local, foreign)
        then:
        transfer1.source == local
        transfer1.target.startsWith(foreign)

        // nothing change then the target path will be the same
        when:
        def transfer2 = porter.createFileTransfer(local, foreign)
        then:
        transfer2.source == local
        transfer2.target == transfer1.target

        when:
        // copy the source to the expected target
        FileHelper.copyPath(local, transfer1.target)
        and:
        // the transfer path should not be modified
        def transfer3 = porter.createFileTransfer(local, foreign)
        then:
        transfer3.source == local
        transfer3.target == transfer1.target

        when:
        // modify the file in the expected target path
        transfer1.target.text = 'Ciao moundo'
        and:
        // the transfer path should not be modified
        def transfer4 = porter.createFileTransfer(local, foreign)
        then:
        transfer4.source == local
        transfer4.target != transfer1.target  // <-- it's changed
        and:
        transfer4.target.startsWith(foreign)  // <-- it's still in the foreign path

        cleanup:
        folder?.deleteDir()
    }
    def 'should stage a file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def local1 = folder.resolve('hola.text')
        def foreign1 = TestHelper.createInMemTempFile('hola.txt', 'hola mundo!')
        and:
        def porter = new FilePorter.FileTransfer(foreign1, local1)

        when:
        porter.stageForeignFile(foreign1, local1)
        then:
        local1.text == foreign1.text

        when:
        def ts = Files.getLastModifiedTime(local1)
        and:
        sleep 1100
        porter.stageForeignFile(foreign1, local1)
        then:
        // file was not touched, since it was in the cache
        ts == Files.getLastModifiedTime(local1)

        cleanup:
        folder?.deleteDir()
    }

    def 'should check valid files' () {
        given:
        def CONTENT = 'hola mundo!'
        def foreign1 = TestHelper.createInMemTempFile('hola.txt', CONTENT)
        and:
        def folder = Files.createTempDirectory('test')
        def local1 = folder.resolve('hola.text'); local1.text = CONTENT
        and:
        def porter = new FilePorter.FileTransfer(foreign1, local1)

        when:
        def equals = porter.checkPathIntegrity(foreign1, local1)
        then:
        equals

        when:
        local1.text = 'foo'
        and:
        equals = porter.checkPathIntegrity(foreign1, local1)
        then:
        !equals

        cleanup:
        folder?.deleteDir()
    }

    def 'should check valid dirs' () {
        given:
        def foreign1 = TestHelper.createInMemTempDir(); foreign1.resolve('sub').mkdir()
        def f1 = foreign1.resolve('file.1'); f1.text = 'file.1'
        def f2 = foreign1.resolve('file.22'); f2.text = 'file.22'
        def f3 = foreign1.resolve('sub/file.333'); f3.text = 'file.333'
        and:
        def local1 = Files.createTempDirectory('test'); local1.resolve('sub').mkdir()
        def l1 = local1.resolve('file.1'); l1.text = 'file.1'
        def l2 = local1.resolve('file.22'); l2.text = 'file.22'
        def l3 = local1.resolve('sub/file.333'); l3.text = 'file.333'
        and:

        when:
        def equals = FilePorter.checkPathIntegrity(foreign1, local1)
        then:
        equals

        when:
        l3.text = 'foo' // <-- change the file size
        and:
        equals = FilePorter.checkPathIntegrity(foreign1, local1)
        then:
        !equals

        cleanup:
        local1?.deleteDir()
    }
}
