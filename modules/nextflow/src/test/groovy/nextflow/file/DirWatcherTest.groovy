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

import static java.nio.file.StandardWatchEventKinds.ENTRY_CREATE
import static java.nio.file.StandardWatchEventKinds.ENTRY_DELETE
import static java.nio.file.StandardWatchEventKinds.ENTRY_MODIFY

import java.nio.file.Files
import java.nio.file.Path

import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DirWatcherTest extends Specification {


    def 'should watch for a new file' () {

        given:
        def folder = Files.createTempDirectory('test')
        def watcher = new DirWatcher('glob', "$folder/", 'hola.txt', false, 'create', folder.getFileSystem())

        when:
        List results = []
        watcher.apply { Path file -> results.add(file.name) }
        sleep 500
        Files.createFile(folder.resolve('hello.txt'))
        Files.createFile(folder.resolve('hola.txt'))
        TestHelper.stopUntil { results.size() == 1  }
        watcher.terminate()

        then:
        results.size() == 1
        results.contains('hola.txt')

        cleanup:
        folder?.deleteDir()

    }


    def 'should watch a path for new files' () {

        given:
        def folder = Files.createTempDirectory('test')
        def watcher = new DirWatcher('glob', "$folder/", '*.txt', false, 'create', folder.getFileSystem())

        when:
        List results = []
        watcher.apply { Path file -> results.add(file.name) }
        sleep 500
        Files.createFile(folder.resolve('hello.txt'))
        Files.createFile(folder.resolve('hola.txt'))
        Files.createFile(folder.resolve('ciao.fasta'))
        TestHelper.stopUntil { results.size() == 2  }
        watcher.terminate()

        then:
        results.size() == 2
        results.contains('hello.txt')
        results.contains('hola.txt')

        cleanup:
        folder?.deleteDir()

    }

    def 'should watch a path for modified and deleted files' () {

        given:
        def folder = Files.createTempDirectory('test')
        Files.createFile(folder.resolve('hello.txt'))
        Files.createFile(folder.resolve('hola.txt'))
        Files.createFile(folder.resolve('ciao.txt'))

        def watcher = new DirWatcher('glob', "$folder/", '*', false, 'modify,delete', folder.getFileSystem())

        when:
        Set results = []
        watcher.apply { Path file -> results.add(file.name) }
        sleep 1000
        folder.resolve('hello.txt').text = '1'
        folder.resolve('ciao.txt').delete()
        TestHelper.stopUntil { results.size() == 2 }
        watcher.terminate()
        
        then:
        results.size() == 2
        results.contains('hello.txt')
        results.contains('ciao.txt')

        cleanup:
        folder?.deleteDir()

    }

    def 'should watch recursively' () {

        given:
        def folder = Files.createTempDirectory('test')
        Files.createDirectories(folder.resolve('foo/bar'))

        def watcher = new DirWatcher('glob', "$folder/", '**.txt', false, 'create', folder.getFileSystem())

        when:
        List results = []
        watcher.apply { Path file -> results.add(folder.relativize(file).toString()) }
        sleep 500
        folder.resolve('hello.txt').text = '1'
        folder.resolve('foo/hola.txt').text = '2'
        folder.resolve('foo/ciao.fa').text = '3'
        folder.resolve('foo/bar/bonjour.txt').text = '4'
        folder.resolve('foo/bar/hallo.fa').text = '5'
        TestHelper.stopUntil { results.size() == 3 }
        watcher.terminate()

        then:
        results.size() == 3
        results.contains('hello.txt')
        results.contains('foo/hola.txt')
        results.contains('foo/bar/bonjour.txt')

        cleanup:
        folder?.deleteDir()

    }

    def 'should watch into subdirs' () {

        given:
        def folder = Files.createTempDirectory('test')
        Files.createDirectories(folder.resolve('foo/bar'))

        def watcher = new DirWatcher('glob', "$folder/", '**/h*', false, 'create', folder.getFileSystem())

        when:
        List results = []
        watcher.apply { Path file -> results.add(folder.relativize(file).toString()) }
        sleep 500
        folder.resolve('hello.txt').text = '1'
        folder.resolve('foo/hola.txt').text = '2'
        folder.resolve('foo/ciao.fa').text = '3'
        folder.resolve('foo/bar/bonjour.txt').text = '4'
        folder.resolve('foo/bar/hallo.fa').text = '5'
        TestHelper.stopUntil { results.size() == 2 }
        watcher.terminate()

        then:
        results.size() == 2
        results.contains('foo/hola.txt')
        results.contains('foo/bar/hallo.fa')

        cleanup:
        folder?.deleteDir()

    }

    def 'should watch into newly create directories' () {

        given:
        def folder = Files.createTempDirectory('test')

        def watcher = new DirWatcher('glob', "$folder/", '**/*.txt', false, 'create', folder.getFileSystem())

        when:
        List results = []
        watcher.apply { Path file -> results.add(folder.relativize(file).toString()) }
        sleep 500
        Files.createDirectories(folder.resolve('foo/bar'))
        sleep 10000
        folder.resolve('hello.txt').text = '1'
        folder.resolve('foo/hola.txt').text = '2'
        folder.resolve('foo/ciao.txt').text = '3'
        folder.resolve('foo/bar/bonjour.txt').text = '4'
        folder.resolve('foo/bar/hallo.doc').text = '5'
        TestHelper.stopUntil { results.size() == 3 }
        watcher.terminate()

        then:
        results.size() == 3
        results.contains('foo/hola.txt')
        results.contains('foo/ciao.txt')
        results.contains('foo/bar/bonjour.txt')

        cleanup:
        folder?.deleteDir()

    }



    def 'should convert string events'() {

        when:
        DirWatcher.stringToWatchEvents('xxx')
        then:
        thrown(IllegalArgumentException)

        expect:
        DirWatcher.stringToWatchEvents() == [ ENTRY_CREATE ]
        DirWatcher.stringToWatchEvents('create,delete') == [ENTRY_CREATE, ENTRY_DELETE]
        DirWatcher.stringToWatchEvents('Create , MODIFY ') == [ENTRY_CREATE, ENTRY_MODIFY]

    }

}
