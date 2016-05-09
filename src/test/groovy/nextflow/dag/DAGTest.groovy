package nextflow.dag

import groovyx.gpars.dataflow.DataflowChannel
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

        def v1=null
        def v2=null

        def dag = new DAG()
        when:
        dag.addVertex(
                DAG.Type.PROCESS,
                'Process 1',
                [ new DAG.ChannelHandler(instance: ch1, label: 'Channel 1') ],
                [ new DAG.ChannelHandler(instance: ch2, label: 'Channel 2') ] )

        v1 = dag.vertices[0]

        then:
        dag.vertices.size() == 1
        v1.label == 'Process 1'
        dag.indexOf(v1) == 0

        dag.edges.size() == 2

        dag.edges[0].label == 'Channel 1'
        dag.edges[0].instance .is ch1
        dag.edges[0].from == null
        dag.edges[0].to == v1

        dag.edges[1].label == 'Channel 2'
        dag.edges[1].instance .is ch2
        dag.edges[1].from == v1
        dag.edges[1].to == null

        when:
        dag.addVertex(
                DAG.Type.PROCESS,
                'Process 2',
                [ new DAG.ChannelHandler(instance: ch2) ],
                [ new DAG.ChannelHandler(instance: ch3, label: 'Channel 3') ] )

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
        dag.edges[0].instance .is ch1
        dag.edges[0].from == null
        dag.edges[0].to == v1

        dag.edges[1].label == 'Channel 2'
        dag.edges[1].instance .is ch2
        dag.edges[1].from == v1
        dag.edges[1].to == v2

        dag.edges[2].label == 'Channel 3'
        dag.edges[2].instance .is ch3
        dag.edges[2].from == v2
        dag.edges[2].to == null


        when:
        dag.addVertex(
                DAG.Type.PROCESS,
                'Process 3',
                [ new DAG.ChannelHandler(instance: ch1, label: 'Channel 1') ],
                null )
        then:
        thrown( MultipleInputChannelException )

        when:
        dag.addVertex(
                DAG.Type.PROCESS,
                'Process 4',
                null,
                [new DAG.ChannelHandler(instance: ch2)] )
        then:
        thrown( MultipleOutputChannelException )
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
                [ new DAG.ChannelHandler(instance: ch1, label: 'Channel 1') ],
                [ new DAG.ChannelHandler(instance: ch2, label: 'Channel 2') ] )

        dag.addVertex(
                DAG.Type.PROCESS,
                'Process 2',
                [ new DAG.ChannelHandler(instance: ch2) ],
                [ new DAG.ChannelHandler(instance: ch3, label: 'Channel 3') ] )

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
                [ new DAG.ChannelHandler(instance: ch1, label: 'Channel 1') ],
                [ new DAG.ChannelHandler(instance: ch2, label: 'Channel 2') ] )


        when:
        dag.resolveEdgeNames(map)

        then:
        dag.edges[0].label == 'channel_1'
        dag.edges[1].label == 'funnel_2'

    }

}
