/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DAGTest extends Specification {

    def 'should create a vertex' () {

        given:
        def dag = new DAG()
        when:
        def v1 = dag.createVertex(DAG.Type.PROCESS, 'Label A')
        def v2 = dag.createVertex(DAG.Type.OPERATOR, 'Label B')

        then:
        v1.label == 'Label A'
        v1.order == 0
        v1.name == 'p0'
        v1.type == DAG.Type.PROCESS

        v2.label == 'Label B'
        v2.order == 1
        v2.name == 'p1'
        v2.type == DAG.Type.OPERATOR
    }


    def 'should add new vertices' () {

        given:
        def ch1 = Mock(DataflowChannel)
        def ch2 = Mock(DataflowChannel)
        def ch3 = Mock(DataflowChannel)

        def v1
        def v2

        def dag = new DAG()
        when:
        dag.addVertex(
                DAG.Type.PROCESS,
                'Process 1',
                [ new DAG.ChannelHandler(channel: ch1, label: 'Channel 1') ],
                [ new DAG.ChannelHandler(channel: ch2, label: 'Channel 2') ] )

        v1 = dag.vertices[0]

        then:
        dag.vertices.size() == 1
        v1.label == 'Process 1'
        dag.indexOf(v1) == 0

        dag.edges.size() == 2

        dag.edges[0].label == 'Channel 1'
        dag.edges[0].channel .is ch1
        dag.edges[0].from == null
        dag.edges[0].to == v1

        dag.edges[1].label == 'Channel 2'
        dag.edges[1].channel .is ch2
        dag.edges[1].from == v1
        dag.edges[1].to == null

        when:
        dag.addVertex(
                DAG.Type.PROCESS,
                'Process 2',
                [ new DAG.ChannelHandler(channel: ch2) ],
                [ new DAG.ChannelHandler(channel: ch3, label: 'Channel 3') ] )

        v1 = dag.vertices[0]
        v2 = dag.vertices[1]
        then:
        dag.vertices.size() == 2
        v1.label == 'Process 1'
        v1.order == 0

        v2.label == 'Process 2'
        v2.order == 1

        dag.edges.size() == 3

        dag.edges[0].label == 'Channel 1'
        dag.edges[0].channel .is ch1
        dag.edges[0].from == null
        dag.edges[0].to == v1

        dag.edges[1].label == 'Channel 2'
        dag.edges[1].channel .is ch2
        dag.edges[1].from == v1
        dag.edges[1].to == v2

        dag.edges[2].label == 'Channel 3'
        dag.edges[2].channel .is ch3
        dag.edges[2].from == v2
        dag.edges[2].to == null

    }

    def 'should throw an exception when the same variable is used multiple times' () {

        given:
        def dag = new DAG()
        def ch1 = new DataflowQueue()
        def ch2 = new DataflowQueue()

        when:
        dag.addVertex( DAG.Type.PROCESS, 'Process 1', [ new DAG.ChannelHandler(channel: ch1) ], null )
        dag.addVertex( DAG.Type.PROCESS, 'Process 2', [ new DAG.ChannelHandler(channel: ch1) ], null )

        then:
        def e1=thrown( MultipleInputChannelException )
        e1.message == 'Channels cannot be used as input in more than one process or operator'

        when:
        dag.addVertex( DAG.Type.PROCESS, 'Process 3', null, [ new DAG.ChannelHandler(channel: ch2) ] )
        dag.addVertex( DAG.Type.PROCESS, 'Process 4', null, [ new DAG.ChannelHandler(channel: ch2) ] )
        then:
        def e2=thrown( MultipleOutputChannelException )
        e2.message == 'Channels cannot be used as output in more than one process or operator'

    }

    def 'should not throw an exception with multiple dataflow input variables' () {
        given:
        def dag = new DAG()
        def ch1 = new DataflowVariable()
        def ch2 = new DataflowVariable()

        when:
        dag.addVertex( DAG.Type.PROCESS, 'Process 1', [ new DAG.ChannelHandler(channel: ch1) ], null )
        dag.addVertex( DAG.Type.PROCESS, 'Process 2', [ new DAG.ChannelHandler(channel: ch1) ], null )
        then:
        dag.vertices.size()==3
        dag.edges.size()==2
        dag.vertices[0].type == DAG.Type.ORIGIN
        // the two edges are the same channel
        dag.edges[0].channel.is ch1
        dag.edges[1].channel.is ch1
        // the two edges share the same origin
        dag.edges[0].from == dag.edges[1].from
        dag.edges[0].from == dag.vertices[0]
        // and end-up to two different vertices
        dag.edges[0].to == dag.vertices[1]
        dag.edges[1].to == dag.vertices[2]

        when:
        dag.addVertex( DAG.Type.PROCESS, 'Process 3', null, [ new DAG.ChannelHandler(channel: ch2) ] )
        dag.addVertex( DAG.Type.PROCESS, 'Process 4', null, [ new DAG.ChannelHandler(channel: ch2) ])
        then:
        def e=thrown(MultipleOutputChannelException)
        e.message == 'Channels cannot be used as output in more than one process or operator'

    }

    def 'should add missing vertices' () {

        given:
        def ch1 = Mock(DataflowChannel)
        def ch2 = Mock(DataflowChannel)
        def ch3 = Mock(DataflowChannel)

        def dag = new DAG()

        when:
        dag.addVertex(
                DAG.Type.PROCESS,
                'Process 1',
                [ new DAG.ChannelHandler(channel: ch1, label: 'Channel 1') ],
                [ new DAG.ChannelHandler(channel: ch2, label: 'Channel 2') ] )

        dag.addVertex(
                DAG.Type.PROCESS,
                'Process 2',
                [ new DAG.ChannelHandler(channel: ch2) ],
                [ new DAG.ChannelHandler(channel: ch3, label: 'Channel 3') ] )

        def p0 = dag.vertices.get(0)
        def p1 = dag.vertices.get(1)
        then:
        dag.vertices.size() == 2
        dag.vertices[0].label == 'Process 1'
        dag.vertices[1].label == 'Process 2'

        dag.edges.size() == 3
        dag.edges[0].from == null
        dag.edges[0].to == p0
        dag.edges[1].from == p0
        dag.edges[1].to == p1
        dag.edges[2].from == p1
        dag.edges[2].to == null

        when:
        dag.normalizeMissingVertices()

        def origin = dag.vertices.get(0)
        def proc1 = dag.vertices.get(1)
        def proc2 = dag.vertices.get(2)
        def term = dag.vertices.get(3)

        then:
        dag.vertices.size() == 4
        dag.vertices[0] == origin
        dag.vertices[0].type == DAG.Type.ORIGIN
        dag.vertices[1].label == 'Process 1'
        dag.vertices[2].label == 'Process 2'
        dag.vertices[3] == term
        dag.vertices[3].type == DAG.Type.NODE

        dag.edges.size() == 3
        dag.edges[0].from == origin
        dag.edges[0].to == proc1
        dag.edges[1].from == proc1
        dag.edges[1].to == proc2
        dag.edges[2].from == proc2
        dag.edges[2].to == term

    }

    def 'should take edge names from variables name' () {
        given:
        def ch1 = Mock(DataflowChannel)
        def ch2 = Mock(DataflowChannel)
        def map = [channel_1: ch1, funnel_2: ch2]

        def dag = new DAG()
        dag.addVertex(
                DAG.Type.PROCESS,
                'Process 1',
                [ new DAG.ChannelHandler(channel: ch1, label: 'Channel 1') ],
                [ new DAG.ChannelHandler(channel: ch2, label: 'Channel 2') ] )


        when:
        dag.resolveEdgeNames(map)

        then:
        dag.edges[0].label == 'channel_1'
        dag.edges[1].label == 'funnel_2'

    }

}
