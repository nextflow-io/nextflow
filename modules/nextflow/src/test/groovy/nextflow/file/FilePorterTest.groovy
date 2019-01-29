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

import spock.lang.Specification

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Paths
import java.util.concurrent.Callable
import java.util.concurrent.ExecutorCompletionService
import java.util.concurrent.Executors

import nextflow.Session
import nextflow.util.Duration
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
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

    def 'should submit actions' () {

        given:
        def folder = Files.createTempDirectory('test')
        def session = new Session(workDir: folder)

        def executor = new ExecutorCompletionService(Executors.newFixedThreadPool(2))
        def actions = [
                { sleep(500) } as Callable,
                { sleep(2000) } as Callable
        ]
        def paths = [Paths.get('/foo/bar')]
        FilePorter porter = Spy(FilePorter, constructorArgs:[session])

        when:
        def futures = porter.submitStagingActions(actions, paths)
        then:
        futures.size() == 2
        0 * porter.getStagingMessage(paths)

        when:
        porter.@pollTimeout = Duration.of('1sec')
        futures = porter.submitStagingActions(actions, paths)
        then:
        futures.size() == 2
        1 * porter.getStagingMessage(paths)
        porter.getStagingMessage(paths) == 'Staging foreign file: /foo/bar'

        cleanup:
        folder?.deleteDir()
    }
}
