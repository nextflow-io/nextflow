import nextflow.processor.NopeTaskProcessor
import nextflow.script.CliRunner
import spock.lang.Specification
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

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FunctionalTest extends Specification {

    /**
     * test passing values through environment variables
     */
    def 'test environment' () {

        setup:
        def environment = [ XXX: 'value1']
        def script = '''

            def x = config.env['XXX']
            def y = config.env['YYY'] ?: -1

            [ x, y ]

            '''

        when:
        def runner = new CliRunner( [env: environment] )

        then:
        runner.execute( script ) == ['value1', -1]

        cleanup:
        runner.workDirectory?.deleteDir()

    }

    /*
     * test passing values through command line parameters
     */
    def 'test parameters' () {

        setup:
        def script = """

            X = 1
            Y = 2
            def Z = 3

            [ X, Y, Z, W ]
            """

        /*
         * the X value is overridden by the 'X' value specified as value
         * the Y value is not overridden because it is not passed as parameter
         * the X value is NOT overridden because is a local variable
         * the W is taken directly from the binding params
         */
        when:
        def runner = new CliRunner()
        runner.setParam( [X:'10', Z: '30', W:'long string']  )
        def result = runner.execute(script)

        then:
        result[0] == 10
        result[1] == 2
        result[2] == 3
        result[3] == 'long string'

        cleanup:
        runner.workDirectory?.deleteDir()

    }

    /*
     * test passing values through command line argument (unnamed parameters)
     */
    def 'test args'()  {

        when:
        def script = """
            def len = args.size()
            def x = args[0]
            def y = args[1]

            return [ x, y, len ]
            """
        def runner = new CliRunner()
        def result = runner.execute(script, 'hello', 'hola' )

        then:
        result[0] == 'hello'
        result[1] == 'hola'
        result[2] == 2

        cleanup:
        runner.workDirectory?.deleteDir()
    }



    def 'test configure processor'() {

        setup:
        def configStr = '''
             task {
                processor = 'nope'
                echo = true
                shell = 'zsh'
                threads = 10
                environment = [a:1, b:2,c:3]
                validExitCodes = [0,11,22,33]
            }
            '''
        def cfg = new ConfigSlurper().parse(configStr)


        when:
        def script = '''

            task('taskHello') { '' }

            return taskProcessor
            '''

        def runner = new CliRunner(cfg)
        def task = runner.execute( script )

        then:
        task instanceof NopeTaskProcessor
        task.getName() == 'taskHello'
        task.getEcho() == true
        task.getShell() == 'zsh'
        task.getThreads() == 10
        task.getEnvironment().entrySet() == [a:1,b:2,c:3].entrySet()
        task.getValidExitCodes() == [0,11,22,33]


        cleanup:
        runner.workDirectory?.deleteDir()

    }



}
