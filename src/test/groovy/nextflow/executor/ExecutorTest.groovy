/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import java.nio.file.Paths

import nextflow.file.FileHolder
import nextflow.processor.ProcessConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ExecutorTest extends Specification {


    def testStagingFilesScript() {
        setup:
        def binding = new Binding()
        def executor = [:] as Executor

        def param1 = new FileInParam(binding, []).bind('file.txt') as FileInParam
        def param2 = new FileInParam(binding, []).bind('seq_*.fa') as FileInParam
        Map<FileInParam,List<FileHolder>> files = [:]
        files[param1] = [FileHolder.get('/home/data/sequences', 'file.txt')]
        files[param2] = [FileHolder.get('/home/data/file1','seq_1.fa'), FileHolder.get('/home/data/file2','seq_2.fa'), FileHolder.get('/home/data/file3','seq_3.fa') ]

        when:
        def script = executor.stagingFilesScript(files)
        def lines = script.readLines()

        then:
        lines[0] == "rm -f 'file.txt'"
        lines[1] == "rm -f 'seq_1.fa'"
        lines[2] == "rm -f 'seq_2.fa'"
        lines[3] == "rm -f 'seq_3.fa'"
        lines[4] == "ln -s '/home/data/sequences' 'file.txt'"
        lines[5] == "ln -s '/home/data/file1' 'seq_1.fa'"
        lines[6] == "ln -s '/home/data/file2' 'seq_2.fa'"
        lines[7] == "ln -s '/home/data/file3' 'seq_3.fa'"
        lines.size() == 8


        when:
        files = [:]
        files[param1] = [FileHolder.get('/data/file', 'file.txt')]
        files[param2] = [FileHolder.get('/data/seq','seq_1.fa') ]
        script = executor.stagingFilesScript(files, '; ')
        lines = script.readLines()

        then:
        lines[0] == "rm -f 'file.txt'; rm -f 'seq_1.fa'; ln -s '/data/file' 'file.txt'; ln -s '/data/seq' 'seq_1.fa'; "
        lines.size() == 1

    }


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

    def testCleanGlobStar() {

        given:
        def executor = [:] as Executor

        expect:
        executor.removeGlobStar('a/b/c') == 'a/b/c'
        executor.removeGlobStar('/a/b/c') == '/a/b/c'
        executor.removeGlobStar('some/*/path') == 'some/*/path'
        executor.removeGlobStar('some/**/path') == 'some'
        executor.removeGlobStar('some/**/path') == 'some'
        executor.removeGlobStar('some*') == 'some*'
        executor.removeGlobStar('some**') == '*'

    }

    def testNormalizePath() {

        given:
        def executor = [:] as Executor

        expect:
        executor.normalizeGlobStarPaths(['file1.txt','path/file2.txt','path/**/file3.txt', 'path/**/file4.txt','**/fa']) == ['file1.txt','path/file2.txt','path','*']
    }

    def 'should return a valid `cp` command' () {

        given:
        def executor = [:] as Executor
        expect:
        executor.copyCommand(source, target) == result

        where:
        source              | target    | result
        'file.txt'          | '/to/dir' | 'cp -fR file.txt /to/dir'
        'path_name'         | '/to/dir' | 'cp -fR path_name /to/dir'
        'input/file.txt'    | '/to/dir' | 'mkdir -p /to/dir/input && cp -fR input/file.txt /to/dir/input'
        'long/path/name'    | '/to/dir' | 'mkdir -p /to/dir/long/path && cp -fR long/path/name /to/dir/long/path'
        'path_name/*'       | '/to/dir' | 'mkdir -p /to/dir/path_name && cp -fR path_name/* /to/dir/path_name'
        'path_name/'        | '/to/dir' | 'mkdir -p /to/dir/path_name && cp -fR path_name/ /to/dir/path_name'

    }

    def 'should return a valid `mv` command' () {

        given:
        def executor = [:] as Executor
        expect:
        executor.copyCommand(source, target, 'move') == result

        where:
        source              | target    | result
        'file.txt'          | '/to/dir' | 'mv -f file.txt /to/dir'
        'file.txt'          | '/to/dir' | 'mv -f file.txt /to/dir'
        'path_name'         | '/to/dir' | 'mv -f path_name /to/dir'
        'input/file.txt'    | '/to/dir' | 'mkdir -p /to/dir/input && mv -f input/file.txt /to/dir/input'
        'long/path/name'    | '/to/dir' | 'mkdir -p /to/dir/long/path && mv -f long/path/name /to/dir/long/path'
        'path_name/*'       | '/to/dir' | 'mkdir -p /to/dir/path_name && mv -f path_name/* /to/dir/path_name'
        'path_name/'        | '/to/dir' | 'mkdir -p /to/dir/path_name && mv -f path_name/ /to/dir/path_name'

    }

    def 'should return a valid `rsync` command' () {

        given:
        def executor = [:] as Executor
        expect:
        executor.copyCommand(source, target, 'rsync') == result

        where:
        source              | target    | result
        'file.txt'          | '/to/dir' | 'rsync -rRl file.txt /to/dir'
        'file.txt'          | '/to/dir' | 'rsync -rRl file.txt /to/dir'
        'path_name'         | '/to/dir' | 'rsync -rRl path_name /to/dir'
        'input/file.txt'    | '/to/dir' | 'rsync -rRl input/file.txt /to/dir'
        'long/path/name'    | '/to/dir' | 'rsync -rRl long/path/name /to/dir'
        'path_name/*'       | '/to/dir' | 'rsync -rRl path_name/* /to/dir'
        'path_name/'        | '/to/dir' | 'rsync -rRl path_name/ /to/dir'

    }

    def 'should return cp script to unstage output files' () {

        given:
        def process = Mock(TaskProcessor)
        process.getConfig() >> Mock(ProcessConfig)

        def executor = [:] as Executor
        def task = Mock(TaskRun)
        task.getProcessor() >> process
        task.getOutputFilesNames() >> [ 'simple.txt', 'my/path/file.bam' ]
        task.getTargetDir() >> Paths.get('/target/work/dir')

        when:
        def script = executor.unstageOutputFilesScript( task )
        then:
        script == '''
                mkdir -p /target/work/dir
                cp -fR simple.txt /target/work/dir || true
                mkdir -p /target/work/dir/my/path && cp -fR my/path/file.bam /target/work/dir/my/path || true
                '''
                .stripIndent().rightTrim()

    }

    def 'should return rsync script to unstage output files' () {

        given:
        def config = new ProcessConfig([unstageStrategy: 'rsync'])

        def process = Mock(TaskProcessor)
        process.getConfig() >> config

        def executor = [:] as Executor
        def task = Mock(TaskRun)
        task.getProcessor() >> process
        task.getOutputFilesNames() >> [ 'simple.txt', 'my/path/file.bam' ]
        task.getTargetDir() >> Paths.get('/target/work/dir')

        when:
        def script = executor.unstageOutputFilesScript( task )
        then:
        script == '''
                mkdir -p /target/work/dir
                rsync -rRl simple.txt /target/work/dir || true
                rsync -rRl my/path/file.bam /target/work/dir || true
                '''
                .stripIndent().rightTrim()

    }

}
