/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.dag

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MultipleOutputChannelExceptionTest extends Specification {

    def 'should return valid error messages' () {

        given:
        def dag = new DAG()
        def node1 = dag.createVertex(DAG.Type.NODE, 'node1')
        def node2 = dag.createVertex(DAG.Type.NODE, 'node2')
        def proc1 = dag.createVertex(DAG.Type.PROCESS, 'x1')
        def proc2 = dag.createVertex(DAG.Type.PROCESS, 'x2')

        expect:
        MultipleOutputChannelException.message(null,null,null) == 'Channels cannot be used as output in more than one process or operator'
        MultipleOutputChannelException.message('channel1',null,null) == 'Channel `channel1` has been used as an output by more than a process or an operator'
        MultipleOutputChannelException.message('channel1',node1,null) == 'Channel `channel1` has been used as an output by more than a process or an operator'
        MultipleOutputChannelException.message('channel1',proc1,null) == 'Channel `channel1` has been used twice as an output by process `x1` and another operator'
        MultipleOutputChannelException.message('channel1',proc1,proc2) == 'Channel `channel1` has been used twice as an output by process `x1` and process `x2`'
        MultipleOutputChannelException.message('channel1',proc1,node2) == 'Channel `channel1` has been used twice as an output by process `x1` and another operator'
        MultipleOutputChannelException.message('channel1',null,node2) == 'Channel `channel1` has been used as an output by more than a process or an operator'
        MultipleOutputChannelException.message('channel1',proc1,proc1) == 'Channel `channel1` has been used twice as an output by process `x1`'
        MultipleOutputChannelException.message('channel1',node1,node1) == 'Channel `channel1` has been used as an output by more than a process or an operator'

    }
    
    
}
