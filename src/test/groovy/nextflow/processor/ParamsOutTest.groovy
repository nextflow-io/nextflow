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
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Nextflow
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ParamsOutTest extends Specification {

    def testOutParam() {

        setup:
        def binding = new Binding()
        def script = Mock(Script)
        script.getBinding() >> binding


        /*
         * creates a new parameter with name 'x'
         * binds to the channel specified by *using*
         */
        when:
        def ch = new DataflowVariable()
        def param = new FileOutParam(script,'x') .to(ch) 
        then:
        param.name == 'x'
        param.outChannel.is(ch)

        /*
         * creates a new *stdout* parameter
         * binds to the channel specified by *using*
         */
        when:
        ch = new DataflowVariable()
        param = new StdOutParam(script) .to(ch) 
        then:
        param.name == '-'
        param.outChannel.is(ch)


        /*
         * test that:
         * - creates a 'FileOutParam'
         * - the channel is specified by its name
         * - since the channel NOT exists in script bindings a new instance is created
         * - the new instance is added to the script bindings
         */
        when:
        def out = new  FileOutParam( script, 'simple.fa' ) .to('channel1') 
        then:
        out.name == 'simple.fa'
        out.outChannel instanceof DataflowQueue
        binding.hasVariable('channel1')
        binding.getVariable('channel1') instanceof DataflowQueue


        /*
         * test that:
         * - creates a 'FileOutParam'
         * - the channel is specified by its name
         * - since the channel NOT exists in script bindings a new instance is created
         * - the new instance is added to the script bindings
         */
        when:
        out = new StdOutParam( script ) .to('channel2') 
        then:
        out.name == '-'
        out.outChannel instanceof DataflowQueue
        binding.hasVariable('channel2')


        /*
         * test that:
         * - the channel is specified by its name
         * - since the channel exists in script bindings a new instance is NOT created
         * - the parameter channel is the SAME instance as the one bound the script bindings
         */
        when:
        def ch3 = Nextflow.val(3)
        binding.setVariable('channel3', ch3)
        out = new FileOutParam(script, 'file.txt') .to('channel3') 
        then:
        out.name == 'file.txt'
        out.outChannel.is(ch3)
        binding.getVariable('channel3').is( ch3 )

        /*
         *
         */
        when:
        def out1 = new FileOutParam( script, 'file.txt' ) .to 'channel4' autoClose false flat true
        def out2 = new FileOutParam( script, 'file.txt' ) .to 'channel4' autoClose true flat false
        def out3 = new ValueOutParam( script, 'x' ) .to 'channel'

        then:
        !(out1 as FileOutParam).autoClose
        (out1 as FileOutParam).flat

        (out2 as FileOutParam).autoClose
        !(out2 as FileOutParam).flat

        out3 instanceof ValueOutParam
        out3.name == 'x'

    }


}
