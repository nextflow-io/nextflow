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

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Session
import nextflow.executor.NopeExecutor
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptInterpreterTest  {

    def 'test createProcessor' () {

        setup:
        def bindings = new Binding([
                input1: 1,
                input2: new DataflowQueue(),
                input3: [1,2,3],
                output1: new DataflowVariable()
                ]
        )
        def session = new Session()
        def nope = new NopeExecutor(session)
        def interpreter = new ScriptInterpreter(bindings,nope)

        when:
        interpreter.defineChannels {
            """
            ${input(input1)}        // 'input1' is a primitive value mapped to a DataflowVariable
            ${input(input2)}        // 'input2' is DataflowQueue
            ${input(input3)}        // 'input3' is a list mapped to a dataflow variable

            ${output(output1)}      // 'output1' is define in the binding
            ${output(output2)}      // 'output2' is NOT declared in the bindings, so it is created on-request
            ${output(output3, input1)}


            ${X.y.z}                // other expression are ingnored
            ${shell HOLA}           // ..
            """
        }


        then:
        // input1 -- was a scalar in the binding, it is stored as a DataflowVariable in the task
        nope.getInput('input1') instanceof DataflowVariable
        nope.getInput('input2') instanceof DataflowQueue
        nope.getInput('input3') instanceof DataflowQueue
        nope.getInput('input4') == null

        nope.getOutput('output1') instanceof DataflowVariable
        nope.getOutput('output2') instanceof DataflowQueue
        nope.getOutput('output3') == null




    }
}
