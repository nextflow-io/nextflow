/*
 * Copyright 2013-2026, Seqera Labs
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

import static java.nio.file.StandardWatchEventKinds.*

import java.nio.file.Files
import java.nio.file.Path

import org.apache.commons.io.monitor.FileAlterationMonitor
import org.apache.commons.io.monitor.FileAlterationObserver
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DirWatcherV2Test extends Specification {

    def 'should watch for a new file' () {

        given:
        def folder = Files.createTempDirectory('test')
        def watcher = newWatcher(folder, 'hola.txt', 'create')

        when:
        List results = []
        watcher.apply { Path file -> results.add(file.name) }
        Files.createFile(folder.resolve('hello.txt'))
        Files.createFile(folder.resolve('hola.txt'))
        TestHelper.stopUntil { results.size() == 1  }

        then:
        results.size() == 1
        results.contains('hola.txt')

        cleanup:
        watcher?.terminate()
        folder?.deleteDir()

    }


    def 'should watch a path for new files' () {

        given:
        def folder = Files.createTempDirectory('test')
        def watcher = newWatcher(folder, '*.txt', 'create')

        when:
        List results = []
        watcher.apply { Path file -> results.add(file.name) }
        Files.createFile(folder.resolve('hello.txt'))
        Files.createFile(folder.resolve('hola.txt'))
        Files.createFile(folder.resolve('ciao.fasta'))
        TestHelper.stopUntil { results.size() == 2  }

        then:
        results.size() == 2
        results.contains('hello.txt')
        results.contains('hola.txt')

        cleanup:
        watcher?.terminate()
        folder?.deleteDir()

    }

    def 'should watch a path for modified and deleted files' () {

        given:
        def folder = Files.createTempDirectory('test')
        Files.createFile(folder.resolve('hello.txt'))
        Files.createFile(folder.resolve('hola.txt'))
        Files.createFile(folder.resolve('ciao.txt'))

        def watcher = newWatcher(folder, '*', 'modify,delete')

        when:
        Set results = []
        watcher.apply { Path file -> results.add(file.name) }
        folder.resolve('hello.txt').text = '1'
        folder.resolve('ciao.txt').delete()
        TestHelper.stopUntil { results.size() == 2 }

        then:
        results.size() == 2
        results.contains('hello.txt')
        results.contains('ciao.txt')

        cleanup:
        watcher?.terminate()
        folder?.deleteDir()

    }

    def 'should watch recursively' () {

        given:
        def folder = Files.createTempDirectory('test')
        Files.createDirectories(folder.resolve('foo/bar'))

        def watcher = newWatcher(folder, '**.txt', 'create')

        when:
        List results = []
        watcher.apply { Path file -> results.add(folder.relativize(file).toString()) }
        folder.resolve('hello.txt').text = '1'
        folder.resolve('foo/hola.txt').text = '2'
        folder.resolve('foo/ciao.fa').text = '3'
        folder.resolve('foo/bar/bonjour.txt').text = '4'
        folder.resolve('foo/bar/hallo.fa').text = '5'
        TestHelper.stopUntil { results.size() == 3 }

        then:
        results.size() == 3
        results.contains('hello.txt')
        results.contains('foo/hola.txt')
        results.contains('foo/bar/bonjour.txt')

        cleanup:
        watcher?.terminate()
        folder?.deleteDir()

    }

    def 'should watch into subdirs' () {

        given:
        def folder = Files.createTempDirectory('test')
        Files.createDirectories(folder.resolve('foo/bar'))

        def watcher = newWatcher(folder, '**/h*', 'create')

        when:
        List results = []
        watcher.apply { Path file -> results.add(folder.relativize(file).toString()) }
        folder.resolve('hello.txt').text = '1'
        folder.resolve('foo/hola.txt').text = '2'
        folder.resolve('foo/ciao.fa').text = '3'
        folder.resolve('foo/bar/bonjour.txt').text = '4'
        folder.resolve('foo/bar/hallo.fa').text = '5'
        TestHelper.stopUntil { results.size() == 2 }

        then:
        results.size() == 2
        results.contains('foo/hola.txt')
        results.contains('foo/bar/hallo.fa')

        cleanup:
        watcher?.terminate()
        folder?.deleteDir()

    }

    def 'should watch into newly create directories' () {

        given:
        def folder = Files.createTempDirectory('test')

        def watcher = newWatcher(folder, '**/*.txt', 'create')

        when:
        List results = []
        watcher.apply { Path file -> results.add(folder.relativize(file).toString()) }
        Files.createDirectories(folder.resolve('foo/bar'))
        folder.resolve('hello.txt').text = '1'
        folder.resolve('foo/hola.txt').text = '2'
        folder.resolve('foo/ciao.txt').text = '3'
        folder.resolve('foo/bar/bonjour.txt').text = '4'
        folder.resolve('foo/bar/hallo.doc').text = '5'
        TestHelper.stopUntil { results.size() == 3 }

        then:
        results.size() == 3
        results.contains('foo/hola.txt')
        results.contains('foo/ciao.txt')
        results.contains('foo/bar/bonjour.txt')

        cleanup:
        watcher?.terminate()
        folder?.deleteDir()

    }



    /** Poll every 100 ms; {@code monitor.start()} snapshots the dir, so no settle sleep is needed. */
    private static DirWatcherV2 newWatcher(Path folder, String pattern, String events) {
        final watcher = new DirWatcherV2('glob', "$folder/", pattern, false, events, folder.getFileSystem())
        watcher.@monitor = new FileAlterationMonitor(100, watcher.@observer as FileAlterationObserver)
        return watcher
    }

    def 'should convert string events'() {
        given:
        def watcher = Spy(DirWatcherV2)
        when:
        watcher.stringToWatchEvents('xxx')
        then:
        thrown(IllegalArgumentException)

        expect:
        watcher.stringToWatchEvents() == [ ENTRY_CREATE ]
        watcher.stringToWatchEvents('create,delete') == [ENTRY_CREATE, ENTRY_DELETE]
        watcher.stringToWatchEvents('Create , MODIFY ') == [ENTRY_CREATE, ENTRY_MODIFY]

    }
}
