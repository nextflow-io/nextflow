/*
 * Copyright (c) 2012, the authors.
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

import java.nio.file.Paths

import nextflow.Session
import nextflow.executor.NopeExecutor
import nextflow.script.BaseScript
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskProcessorTest extends Specification {


    def testEnvironment() {

        setup:
        def home = File.createTempDir()
        def binFolder = new File(home,'bin')
        binFolder.mkdirs()

        when:
        def wrapper = new DummyScript()
        def session = new Session([env: [X:"1", Y:"2"]])
        session.setBaseDir(home)
        def processor = new DummyProcessor(session, wrapper, new TaskConfig(wrapper))
        def builder = new ProcessBuilder()
        builder.environment().putAll( processor.getProcessEnvironment() )

        then:
        noExceptionThrown()
        builder.environment().X == '1'
        builder.environment().Y == '2'
        builder.environment().PATH == binFolder.toString()

        when:
        session = new Session([env: [X:"1", Y:"2", PATH:'/some']])
        session.setBaseDir(home)
        processor = new DummyProcessor(session, wrapper, new TaskConfig(wrapper))
        builder = new ProcessBuilder()
        builder.environment().putAll( processor.getProcessEnvironment() )

        then:
        noExceptionThrown()
        builder.environment().X == '1'
        builder.environment().Y == '2'
        builder.environment().PATH == "${binFolder.toString()}:/some"


        cleanup:
        home.deleteDir()

    }

    def testFetchInterpreter() {

        when:
        def wrapper = new DummyScript()
        def processor = new DummyProcessor(new Session(), wrapper, new TaskConfig(wrapper))
        def script =
            '''
            #!/bin/perl
            do this
            do that
            '''
        def i = processor.fetchInterpreter(script.stripIndent().trim())
        then:
        i == '/bin/perl'

        when:
        processor = new DummyProcessor(new Session(), wrapper, new TaskConfig(wrapper))
        i = processor.fetchInterpreter('do this')
        then:
        i == null
    }


    def testWildcardExpand() {

        setup:
        def wrapper = new DummyScript()
        def processor = new DummyProcessor(new Session(), wrapper, new TaskConfig(wrapper))

        /*
         * The name do not contain any wildcards *BUT* when multiple files are provide
         * an index number is added to the specified name
         */
        when:
        def list1 = processor.expandWildcards('file_name', Paths.get('x'))
        def list2 = processor.expandWildcards('file_name', [Paths.get('x'), Paths.get('y')])
        then:
        list1 == ['file_name']
        list2 == ['file_name1', 'file_name2']


        /*
         * The star wildcard: when a single item is provided, it is simply ignored
         * When a collection of files is provided, the name is expanded to the index number
         */
        when:
        list1 = processor.expandWildcards('file*.fa', Paths.get('x'))
        list2 = processor.expandWildcards('file_*.fa', [Paths.get('x'), Paths.get('y'), Paths.get('z')])
        then:
        list1 == ['file.fa']
        list2 == ['file_1.fa', 'file_2.fa', 'file_3.fa']

        /*
         * The question mark wildcards *always* expand to an index number
         */
        when:
        list1 = processor.expandWildcards('file?.fa', Paths.get('0'))
        list2 = processor.expandWildcards('file_???.fa', 1..4 )
        def list3 = processor.expandWildcards('file_?.fa', 1..12 )
        then:
        list1 == ['file1.fa']
        list2 == ['file_001.fa', 'file_002.fa', 'file_003.fa', 'file_004.fa']
        list3 == ['file_1.fa', 'file_2.fa', 'file_3.fa', 'file_4.fa', 'file_5.fa', 'file_6.fa', 'file_7.fa', 'file_8.fa', 'file_9.fa', 'file_10.fa', 'file_11.fa', 'file_12.fa']

        when:
        list1 = processor.expandWildcards('*', Paths.get('a'))
        list2 = processor.expandWildcards('*', [Paths.get('x'), Paths.get('y'), Paths.get('z')])
        then:
        list1 == ['a']
        list2 == ['x','y','z']

        when:
        processor.expandWildcards('*', 0)
        then:
        thrown(IllegalArgumentException)

    }


    def testStagingFilesScript() {
        setup:
        def owner = Mock(Script)
        def wrapper = new DummyScript()
        def processor = new DummyProcessor(new Session(), wrapper, new TaskConfig(wrapper))
        def param1 = new FileInParam(owner, 'file.txt')
        def param2 = new FileInParam(owner, 'seq_*.fa')
        Map<FileInParam,Object> files = [:]
        files[param1] = Paths.get('/home/data/sequences')
        files[param2] = [Paths.get('/home/data/file1'), Paths.get('/home/data/file2'), Paths.get('/home/data/file3') ]

        when:
        def script = processor.stagingFilesScript(files)
        def lines = script.readLines()

        then:
        lines[0] == 'rm -f file.txt'
        lines[1] == 'rm -f seq_1.fa'
        lines[2] == 'rm -f seq_2.fa'
        lines[3] == 'rm -f seq_3.fa'
        lines[4] == 'ln -s /home/data/sequences file.txt'
        lines[5] == 'ln -s /home/data/file1 seq_1.fa'
        lines[6] == 'ln -s /home/data/file2 seq_2.fa'
        lines[7] == 'ln -s /home/data/file3 seq_3.fa'
        lines.size() == 8

    }



    static class DummyProcessor extends TaskProcessor {

        DummyProcessor(Session session, BaseScript script, TaskConfig taskConfig) {
            super(new NopeExecutor(), session, new DummyScript(), taskConfig, {})
        }

        @Override
        protected void createOperator() {
        }
    }

    static class DummyScript extends BaseScript {
        @Override
        Object run() {
            return null
        }
    }

}
