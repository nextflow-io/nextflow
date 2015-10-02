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

import java.nio.file.Paths

import nextflow.processor.ProcessConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SimpleFileCopyStrategyTest extends Specification {


    def 'should return stage input file script'() {
        setup:
        def files = [:]
        files['file.txt'] = Paths.get('/home/data/sequences')
        files['seq_1.fa'] = Paths.get('/home/data/file1')
        files['seq_2.fa'] = Paths.get('/home/data/file2')
        files['seq_3.fa'] = Paths.get('/home/data/file3')

        when:
        def strategy = new SimpleFileCopyStrategy(inputFiles: files)
        def script = strategy.getStageInputFilesScript()
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
        files['file.txt'] = Paths.get('/data/file')
        files['seq_1.fa'] = Paths.get('/data/seq')
        strategy = new SimpleFileCopyStrategy(inputFiles: files, separatorChar: '; ')
        script = strategy.getStageInputFilesScript()
        lines = script.readLines()

        then:
        lines[0] == "rm -f 'file.txt'; rm -f 'seq_1.fa'; ln -s '/data/file' 'file.txt'; ln -s '/data/seq' 'seq_1.fa'; "
        lines.size() == 1

    }


    def 'should remove star glob pattern'() {

        given:
        def strategy = [:] as SimpleFileCopyStrategy

        expect:
        strategy.removeGlobStar('a/b/c') == 'a/b/c'
        strategy.removeGlobStar('/a/b/c') == '/a/b/c'
        strategy.removeGlobStar('some/*/path') == 'some/*/path'
        strategy.removeGlobStar('some/**/path') == 'some'
        strategy.removeGlobStar('some/**/path') == 'some'
        strategy.removeGlobStar('some*') == 'some*'
        strategy.removeGlobStar('some**') == '*'

    }

    def 'should normalize path'() {

        given:
        def strategy = [:] as SimpleFileCopyStrategy

        expect:
        strategy.normalizeGlobStarPaths(['file1.txt','path/file2.txt','path/**/file3.txt', 'path/**/file4.txt','**/fa']) == ['file1.txt','path/file2.txt','path','*']
    }

    def 'should return a valid `cp` command' () {

        given:
        def strategy = [:] as SimpleFileCopyStrategy
        expect:
        strategy.copyCommand(source, target) == result

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
        def strategy = [:] as SimpleFileCopyStrategy
        expect:
        strategy.copyCommand(source, target, 'move') == result

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
        def strategy = [:] as SimpleFileCopyStrategy
        expect:
        strategy.copyCommand(source, target, 'rsync') == result

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

        def task = Mock(TaskRun)
        task.getProcessor() >> process
        task.getOutputFilesNames() >> [ 'simple.txt', 'my/path/file.bam' ]
        task.getTargetDir() >> Paths.get('/target/work/dir')

        when:
        def strategy = new SimpleFileCopyStrategy(task)
        def script = strategy.getUnstageOutputFilesScript()
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

        def task = Mock(TaskRun)
        task.getProcessor() >> process
        task.getOutputFilesNames() >> [ 'simple.txt', 'my/path/file.bam' ]
        task.getTargetDir() >> Paths.get('/target/work/dir')

        when:
        def strategy = new SimpleFileCopyStrategy(task)
        def script = strategy.getUnstageOutputFilesScript()
        then:
        script == '''
                mkdir -p /target/work/dir
                rsync -rRl simple.txt /target/work/dir || true
                rsync -rRl my/path/file.bam /target/work/dir || true
                '''
                .stripIndent().rightTrim()

    }

}
