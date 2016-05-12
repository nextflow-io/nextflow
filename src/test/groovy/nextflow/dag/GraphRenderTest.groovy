package nextflow.dag

import groovyx.gpars.dataflow.DataflowChannel
import spock.lang.Specification
import java.nio.file.Files
import nextflow.Session
import nextflow.script.InputsList
import nextflow.script.OutputsList

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Mike Smoot <mes@aescon.com>
 */
class GraphRenderTest extends Specification {

    def DAG test_dag

    def setup() {

        def src = Mock(DataflowChannel)
        def ch1 = Mock(DataflowChannel)
        def ch2 = Mock(DataflowChannel)
        def ch3 = Mock(DataflowChannel)

        test_dag = new DAG()

        test_dag.addVertex(
                DAG.Type.ORIGIN,
                'Source',
                null,
                [ new DAG.ChannelHandler(instance: src, label: 'src') ] )

        test_dag.addVertex(
                DAG.Type.PROCESS,
                'Process 1',
                [ new DAG.ChannelHandler(instance: src, label: 'Source') ],
                [ new DAG.ChannelHandler(instance: ch1, label: 'Channel 1') ] )

        test_dag.addVertex(
                DAG.Type.OPERATOR,
                'Filter',
                [ new DAG.ChannelHandler(instance: ch1, label: 'Channel 1') ],
                [ new DAG.ChannelHandler(instance: ch2, label: 'Channel 2') ] )


        test_dag.addVertex(
                DAG.Type.PROCESS,
                'Process 2',
                [ new DAG.ChannelHandler(instance: ch2, label: 'Channel 2') ],
                [ new DAG.ChannelHandler(instance: ch3, label: 'Channel 3') ] )

    }

    def 'should write a dot file' () {
        given:
        def dot = Files.createTempFile('nxf-','.dot')
        def gr = new GraphRender(dot)
        gr.dag = test_dag

        when:
        gr.onFlowComplete()

        then:
        def result = dot.text
        // is dot
        result.contains('digraph graphname {')

        // contains expected nodes
        result.contains('label="Source"')
        result.contains('label="Process 1"')
        result.contains('label="Filter"')
        result.contains('label="Process 2"')

        // contains at least one edge
        result.contains('->')
    }

    def 'should write a json file' () {
        given:
        def json = Files.createTempFile('nxf-','.json')
        def gr = new GraphRender(json)
        gr.dag = test_dag

        when:
        gr.onFlowComplete()

        then:
        def result = json.text

        // is json
        result.contains('elements: {')
        result.contains('nodes: [')
        result.contains('edges: [')

        // contains expected nodes
        result.contains("label: 'Source'")
        result.contains("label: 'Process 1'")
        result.contains("label: 'Filter'")
        result.contains("label: 'Process 2'")

        result.contains("classes: 'ORIGIN'")
        result.contains("classes: 'PROCESS'")
        result.contains("classes: 'OPERATOR'")

        // contains at least one edge
        result.contains('source:')
        result.contains('target:')
    }

    def 'should write an html file' () {
        given:
        def html = Files.createTempFile('nxf-','.html')
        def gr = new GraphRender(html)
        gr.dag = test_dag

        when:
        gr.onFlowComplete()

        then:
        def result = html.text

        // is html
        result.contains('<html>')
        result.contains('</html>')

        // contains some of the expected json
        result.contains("label: 'Source'")
        result.contains("label: 'Process 1'")
        result.contains("label: 'Filter'")
        result.contains("label: 'Process 2'")
    }

    def 'should write an htm file' () {
        given:
        def htm = Files.createTempFile('nxf-','.htm')
        def gr = new GraphRender(htm)
        gr.dag = test_dag

        when:
        gr.onFlowComplete()

        then:
        def result = htm.text

        // is htm
        result.contains('<html>')
        result.contains('</html>')
        result.contains('<svg')

        // contains some of the expected javascript
        result.contains("label: 'Source'")
        result.contains("label: 'Process 1'")
        result.contains("label: 'Filter'")
        result.contains("label: 'Process 2'")
    }

    def 'should write an svg file' () {
        given:
        def svg = Files.createTempFile('nxf-','.svg')
        def gr = new GraphRender(svg)
        gr.dag = test_dag

        when:
        gr.onFlowComplete()

        then:
        // just assert that _something_ has been written
        svg.size() > 0
    }

    def 'should write a png file' () {
        given:
        def png = Files.createTempFile('nxf-','.png')
        def gr = new GraphRender(png)
        gr.dag = test_dag

        when:
        gr.onFlowComplete()

        then:
        // just assert that _something_ has been written
        png.size() > 0
    }

    def 'should write a pdf file' () {
        given:
        def pdf = Files.createTempFile('nxf-','.pdf')
        def gr = new GraphRender(pdf)
        gr.dag = test_dag

        when:
        gr.onFlowComplete()

        then:
        // just assert that _something_ has been written
        pdf.size() > 0
    }

    def 'should not write anything with unknown suffix' () {
        given:
        def nope = Files.createTempFile('nxf-','.nope')
        def gr = new GraphRender(nope)
        gr.dag = test_dag

        when:
        gr.onFlowComplete()

        then:
        // nothing should be written!
        nope.size() == 0

        // Ideally also assert that an error message has been
        // written to the logs.
    }
}
