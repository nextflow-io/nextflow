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

import nextflow.script.FileInParam
import nextflow.processor.FileHolder
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


}
