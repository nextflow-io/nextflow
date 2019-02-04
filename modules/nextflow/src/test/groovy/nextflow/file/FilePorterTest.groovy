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

import spock.lang.Ignore
import spock.lang.Specification

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.Callable
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit
import java.util.concurrent.TimeoutException

import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.ProcessStageException
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class FilePorterTest extends Specification {

    def 'should get the core threads value' () {

        when:
        def session = new Session()
        def porter = new FilePorter(session)
        then:
        porter.coreThreads == Runtime.getRuntime().availableProcessors()

        when:
        session = new Session([filePorter: [coreThreads: 10]])
        porter = new FilePorter(session)
        then:
        porter.coreThreads == 10

    }

    def 'should get the max threads value' () {

        when:
        def session = new Session()
        def porter = new FilePorter(session)
        then:
        porter.maxThreads == 2 * Runtime.getRuntime().availableProcessors()

        when:
        session = new Session([filePorter: [maxThreads: 99]])
        porter = new FilePorter(session)
        then:
        porter.maxThreads == 99

    }

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
        def local = Paths.get('local.txt')
        def files = [foo: local, bar: foreign1, baz: foreign2]

        when:
        def porter = new FilePorter(session)
        def result = porter.stageForeignFiles(files, folder)
        then:
        result.foo ==  Paths.get('local.txt')

        result.bar.name == 'hola.txt'
        result.bar.text == 'hola mundo!'
        result.bar.fileSystem == FileSystems.default

        result.baz.name == 'ciao.txt'
        result.baz.text == 'ciao mondo!'
        result.baz.fileSystem == FileSystems.default

        cleanup:
        folder?.deleteDir()
    }


    static class DummyStage extends FilePorter.FileStageAction {
        int secs
        DummyStage(String name, Path file, int secs ) {
            super(name,file,Paths.get('/work/dir'),0)
            this.secs = secs
        }

        @Override
        FilePorter.NamePathPair call() throws Exception {
            log.debug "Dummy staging $path"
            sleep secs
            return new FilePorter.NamePathPair(name, path.resolve(name))
        }
    }

    static class ErrorStage extends FilePorter.FileStageAction {

        @Override
        FilePorter.NamePathPair call() throws Exception {
            throw new ProcessStageException('Cannot stage gile')
        }
    }

    def 'should submit actions' () {

        given:
        def folder = Files.createTempDirectory('test')
        def session = new Session(workDir: folder)
        def porter = new FilePorter(session)

        when:
        def actions = [
                new DummyStage('foo', Paths.get('/data/foo'), 500),
                new DummyStage('bar', Paths.get('/data/bat'), 3000)
        ]

        def futures = porter.submitStagingActions(actions)
        then:
        futures.size() == 2
        futures[0].isDone()
        futures[1].isDone()

        when:
        porter.submitStagingActions([ new ErrorStage() ])
        then:
        thrown(ProcessStageException)

        cleanup:
        folder?.deleteDir()
    }


    @Ignore
    def 'should access future' () {

        given:
        def exec = Executors.newFixedThreadPool(2)
        when:
        def fut = exec.submit( { log.info 'start'; sleep 5000; log.info 'done'; return 'Hello' } as Callable )

        def wait = {
            while( true ) {
                try {
                    def str = fut.get(1, TimeUnit.SECONDS)
                    log.info "message => $str"
                    break
                }
                catch( TimeoutException e ) {
                    //
                    log.info('timeout')
                }
            }
        }

        def t1 = Thread.start(wait)
        def t2 = Thread.start(wait)

        t1.join()
        t2.join()

        then:
        noExceptionThrown()
    }
}
