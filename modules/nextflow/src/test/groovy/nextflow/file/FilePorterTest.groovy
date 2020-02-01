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

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.ProcessStageException
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
}
