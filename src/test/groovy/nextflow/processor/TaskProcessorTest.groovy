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

import nextflow.Session
import spock.lang.Specification
import test.DummyProcessor
import test.DummyScript

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



}
