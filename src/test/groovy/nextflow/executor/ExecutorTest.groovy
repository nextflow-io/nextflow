/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

package nextflow.executor

import java.nio.file.Files

import nextflow.script.FileOutParam
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ExecutorTest extends Specification {


    def testCollectOutputFiles() {

        given:
        def param
        def result
        def executor = [:] as Executor

        def folder = Files.createTempDirectory('test')
        folder.resolve('file1.txt').text = 'file 1'
        folder.resolve('file2.fa').text = 'file 2'
        folder.resolve('.hidden.fa').text = 'hidden'
        folder.resolve('dir1').mkdir()
        folder.resolve('dir1').resolve('file3.txt').text = 'file 3'
        folder.resolve('dir1')
        folder.resolve('dir1').resolve('dir2').mkdirs()
        folder.resolve('dir1').resolve('dir2').resolve('file4.fa').text = 'file '
        Files.createSymbolicLink( folder.resolve('dir_link'), folder.resolve('dir1') )

        when:
        result = executor.collectResultFile(folder, '*.fa', 'test', Mock(FileOutParam) )
        then:
        result.collect { it.name }  == ['file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.type('file')
        result = executor.collectResultFile(folder, '*.fa', 'test', param)
        then:
        result.collect { it.name }  == ['file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.type('dir')
        result = executor.collectResultFile(folder, '*.fa', 'test', param)
        then:
        result == []

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        result = executor.collectResultFile(folder, '**.fa', 'test', param)
        then:
        result.collect { it.name }.sort()  == ['file2.fa','file4.fa','file4.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.followLinks(false)
        result = executor.collectResultFile(folder, '**.fa', 'test', param)
        then:
        result.collect { it.name }.sort()  == ['file2.fa','file4.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.maxDepth(1)
        result = executor.collectResultFile(folder, '**.fa', 'test', param)
        then:
        result.collect { it.name }.sort()  == ['file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        result = executor.collectResultFile(folder, '*', 'test', param)
        then:
        result.collect { it.name }.sort()  == ['dir1', 'dir_link', 'file1.txt', 'file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.type('dir')
        result = executor.collectResultFile(folder, '*', 'test', param)
        then:
        result.collect { it.name }.sort()  == ['dir1', 'dir_link']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.type('file')
        result = executor.collectResultFile(folder, '*', 'test', param)
        then:
        result.collect { it.name }.sort()  == ['file1.txt', 'file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.type('file')
        param.hidden(true)
        result = executor.collectResultFile(folder, '*', 'test', param)
        then:
        result.collect { it.name }.sort()  == ['.hidden.fa', 'file1.txt', 'file2.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        result = executor.collectResultFile(folder,'.*', 'test', param)
        then:
        result.collect { it.name }.sort()  == ['.hidden.fa']

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        result = executor.collectResultFile(folder,'file{1,2}.{txt,fa}', 'test', param)
        then:
        result.collect { it.name }.sort() == ['file1.txt', 'file2.fa']

        cleanup:
        folder?.deleteDir()

    }

    def testCollectResultOpts() {

        given:
        def param
        def executor = [:] as Executor

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        then:
        executor.collectResultOpts(param,'file.txt') == [type:'any', followLinks: true, maxDepth: null, hidden: false, relative: false]
        executor.collectResultOpts(param,'path/**') == [type:'file', followLinks: true, maxDepth: null, hidden: false, relative: false]
        executor.collectResultOpts(param,'.hidden_file') == [type:'any', followLinks: true, maxDepth: null, hidden: true, relative: false]

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.type('dir')
        then:
        executor.collectResultOpts(param,'dir-name') == [type:'dir', followLinks: true, maxDepth: null, hidden: false, relative: false]

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.hidden(true)
        then:
        executor.collectResultOpts(param,'dir-name') == [type:'any', followLinks: true, maxDepth: null, hidden: true, relative: false]

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.followLinks(false)
        then:
        executor.collectResultOpts(param,'dir-name') == [type:'any', followLinks: false, maxDepth: null, hidden: false, relative: false]

        when:
        param = new FileOutParam(Mock(Binding), Mock(List))
        param.maxDepth(5)
        then:
        executor.collectResultOpts(param,'dir-name') == [type:'any', followLinks: true, maxDepth: 5, hidden: false, relative: false]
    }

}
