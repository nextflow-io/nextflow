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

package nextflow.processor

import spock.lang.Ignore
import spock.lang.Specification

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import nextflow.util.KryoHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskPathTest extends Specification {

    def 'should create a task path' () {

        given:
        def folder = Files.createTempDirectory('test')
        def hello = folder.resolve('hello.txt')
        hello.text = 'Hello world'
        
        when:
        def path = new TaskPath(hello, 'ciao.txt')
        then:
        path instanceof Path
        "$path" == 'ciao.txt'
        path.getName() == 'ciao.txt'
        path.getFileName() == Paths.get('ciao.txt')
        path.getSimpleName() == 'ciao'
        path.getBaseName() == 'ciao'
        path.getExtension() == 'txt'
        path.getNameCount() == 1
        path.getName(0) == Paths.get('ciao.txt')
        path.startsWith('ciao') == Paths.get('ciao.txt').startsWith('ciao')
        path.startsWith('ciao.txt') == Paths.get('ciao.txt').startsWith('ciao.txt')
        path.normalize() == Paths.get('ciao.txt')
        path.resolve('foo.bar') == Paths.get('ciao.txt').resolve('foo.bar')
        path.resolveSibling('foo.bar') == Paths.get('foo.bar')
        !path.empty()
        path.exists()
        path.size() == hello.size()
        path.lastModified() == hello.lastModified()
        !path.isDirectory()
        path.isFile()
        path.isLink()
        path.canRead()
        path.canWrite()
        !path.canExecute()
        path.permissions == hello.permissions
        path.resolveSymLink() == hello
        path.toRealPath() == hello
        path.toFile() == new File('ciao.txt')
        path.iterator().next() == path.getFileName()
        // ? path.readAttributes() == FilesEx.readAttributes(hello)
        path.matches('ciao*')
        path.getUri() == Paths.get('ciao.txt').toUri()

        // some math
        path + '_bak' == Paths.get('ciao.txt')  + '_bak'
        path / 'foo'  == Paths.get('ciao.txt') / 'foo'
        path || 'foo' == Paths.get('ciao.txt') || 'foo'
        path - 1 == null

        !path.isAbsolute()
        path.getRoot() == null
        path.getParent() == null
        path.getFileSystem() == FileSystems.default

        cleanup:
        folder.deleteDir()

    }

    def 'should validate equals method' () {
        given:
        def p1 = Paths.get('/foo/bar/baz.txt')
        def p2 = Paths.get('/xxx/yyy/zzz.txt')
        def t1 = new TaskPath(p1)
        def t2 = new TaskPath(p1, 'foo.txt')
        def t3 = new TaskPath(p2)

        expect: 
        t1.equals(t1)
        t1.equals(new TaskPath(p1))
        new TaskPath(p1).equals(t1)

        !t1.equals(t2)
        !t2.equals(t1)

        t2.equals(t2)
        t2.equals(new TaskPath(p1, 'foo.txt'))
        new TaskPath(p1, 'foo.txt').equals(t2)

        !t1.equals(t3)

    }


    def 'should validate operator equality' () {
        given:
        def p1 = Paths.get('/foo/bar/baz.txt')
        def p2 = Paths.get('/xxx/yyy/zzz.txt')
        def t1 = new TaskPath(p1)
        def t2 = new TaskPath(p1, 'foo.txt')
        def t3 = new TaskPath(p2)
        def t4 = new TaskPath(p2)
        
        expect:
        TaskPath.equals(p1, t1)
        TaskPath.equals(t1, p1)

        !TaskPath.equals(p1, p2)
        !TaskPath.equals(t1, t2)

        TaskPath.equals(t3,t4)
    }

    def 'should return size zero for non-existing file' () {

        expect:
        new TaskPath(Paths.get('/foo')).size() == 0

    }

    @Ignore
    def 'should report error' () {

        given:
        def folder = Files.createTempDirectory('test')
        def hello = folder.resolve('hello.txt')
        hello.text = 'Hello world'

        when:
        def path = new TaskPath(hello, 'ciao.txt')

        then:
        // -- PROBLEMS
        path == hello  // it fails because Path implements Comparable https://stackoverflow.com/questions/28355773/in-groovy-why-does-the-behaviour-of-change-for-interfaces-extending-compar#comment45123447_28387391
        path.text == 'Hello world'

        cleanup:
        folder.deleteDir()

    }


    def 'should serialised task path' () {

        given:
        def p = new TaskPath(Paths.get('/foo/bar/bax.txt'), 'custom-name.fasta')
        when:
        def buffer = KryoHelper.serialize(p)
        def copy = (TaskPath)KryoHelper.deserialize(buffer)
        then:
        copy.equals(p)
    }

}



