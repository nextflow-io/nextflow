/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

package nextflow.dag

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MultipleInputChannelExceptionTest extends Specification {

    def 'should return valid error messages' () {

        given:
        def dag = new DAG()
        def node1 = dag.createVertex(DAG.Type.NODE, 'node1')
        def node2 = dag.createVertex(DAG.Type.NODE, 'node2')
        def proc1 = dag.createVertex(DAG.Type.PROCESS, 'x1')
        def proc2 = dag.createVertex(DAG.Type.PROCESS, 'x2')

        expect:
        MultipleInputChannelException.message(null,null,null) == 'Channels cannot be used as input in more than one process or operator'
        MultipleInputChannelException.message('channel1',null,null) == 'Channel `channel1` has been used as an input by more than a process or an operator'
        MultipleInputChannelException.message('channel1',node1,null) == 'Channel `channel1` has been used as an input by more than a process or an operator'
        MultipleInputChannelException.message('channel1',proc1,null) == 'Channel `channel1` has been used twice as an input by process `x1` and another operator'
        MultipleInputChannelException.message('channel1',proc1,proc2) == 'Channel `channel1` has been used twice as an input by process `x1` and process `x2`'
        MultipleInputChannelException.message('channel1',proc1,node2) == 'Channel `channel1` has been used twice as an input by process `x1` and another operator'
        MultipleInputChannelException.message('channel1',null,node2) == 'Channel `channel1` has been used as an input by more than a process or an operator'
        MultipleInputChannelException.message('channel1',proc1,proc1) == 'Channel `channel1` has been used twice as an input by process `x1`'
        MultipleInputChannelException.message('channel1',node1,node1) == 'Channel `channel1` has been used as an input by more than a process or an operator'

    }

}
