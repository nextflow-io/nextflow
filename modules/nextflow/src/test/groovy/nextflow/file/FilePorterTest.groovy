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

import nextflow.Global
import nextflow.Session
import test.TestHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FilePorterTest extends Specification {

    def 'should get the core threads value' () {

        when:
        new Session()
        then:
        FilePorter.getCoreThreads() == Runtime.getRuntime().availableProcessors()

        when:
        new Session([filePorter: [coreThreads: 10]])
        then:
        FilePorter.getCoreThreads() == 10

    }

    def 'should get the max threads value' () {

        when:
        new Session()
        then:
        FilePorter.getMaxThreads() == 2 * Runtime.getRuntime().availableProcessors()

        when:
        new Session([filePorter: [maxThreads: 99]])
        then:
        FilePorter.getMaxThreads() == 99

    }

    def 'should get the max retries value' () {

        when:
        new Session()
        then:
        FilePorter.getMaxRetries() == 3

        when:
        new Session([filePorter: [maxRetries: 88]])
        then:
        FilePorter.getMaxRetries() == 88

    }

    def 'should get the download path' () {

        given:
        new Session()
        def fp = Spy(FilePorter)

        when:
        fp.stagingDir
        then:
        1 * fp.getStagingDir() >> Global.session.workDir.resolve('download')

        when:
        Global.session.config = [filePorter: [downloadPath: '/foo/bar']]
        fp.stagingDir
        then:
        1 * fp.getStagingDir() >> Paths.get('/foo/bar')

    }

    def 'should copy foreign files' () {

        given:
        new Session()
        def folder = Files.createTempDirectory('test')
        def foreign1 = TestHelper.createInMemTempFile('hola.txt', 'hola mundo!')
        def foreign2 = TestHelper.createInMemTempFile('ciao.txt', 'ciao mondo!')
        def local = Paths.get('local.txt')
        def files = [foo: local, bar: foreign1, baz: foreign2]

        when:
        def d = new FilePorter(folder)
        def result = d.stageForeignFiles(files)
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
        new Session()
        def executor = new ExecutorCompletionService(Executors.newFixedThreadPool(2))
        def actions = [
                { sleep(500) } as Callable,
                { sleep(2000) } as Callable
        ]
        def paths = [Paths.get('/foo/bar')]
        def fp = Spy(FilePorter)

        when:
        def futures = fp.submitStagingActions(executor, actions, paths)
        then:
        futures.size() == 2
        0 * fp.getDownloadMessage(paths)

        when:
        fp.pollTimeout = 1L
        futures = fp.submitStagingActions(executor, actions, paths)
        then:
        futures.size() == 2
        1 * fp.getDownloadMessage(paths)
        fp.getDownloadMessage(paths) == 'Staging foreign file /foo/bar'
    }
}
