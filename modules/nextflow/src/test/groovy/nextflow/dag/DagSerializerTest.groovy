/*
 * Copyright 2013-2026, Seqera Labs
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
import nextflow.Session
import spock.lang.Specification
/**
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
class DagSerializerTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def 'should return null for a null dag' () {
        expect:
        DagSerializer.toMap(null) == null
    }

    def 'should serialize an empty dag' () {
        when:
        def result = DagSerializer.toMap(new DAG())
        then:
        result == [vertices: [], edges: []]
    }

    def 'should faithfully serialize vertices and edges' () {
        given:
        def ch1 = Mock(DataflowChannel)
        def ch2 = Mock(DataflowChannel)
        def ch3 = Mock(DataflowChannel)
        and:
        def dag = new DAG()
        // PROCESS --ch2--> OPERATOR --ch3--> PROCESS
        dag.addVertex(
                DAG.Type.PROCESS,
                'foo',
                [ new DAG.ChannelHandler(channel: ch1, label: 'in') ],
                [ new DAG.ChannelHandler(channel: ch2, label: 'out') ] )
        dag.addVertex(
                DAG.Type.OPERATOR,
                'map',
                [ new DAG.ChannelHandler(channel: ch2) ],
                [ new DAG.ChannelHandler(channel: ch3) ] )
        dag.addVertex(
                DAG.Type.PROCESS,
                'bar',
                [ new DAG.ChannelHandler(channel: ch3) ],
                [] )

        when:
        def result = DagSerializer.toMap(dag)
        def vFoo = result.vertices.find { it.label == 'foo' }
        def vMap = result.vertices.find { it.label == 'map' }
        def vBar = result.vertices.find { it.label == 'bar' }

        then:
        // every declared vertex is preserved, including the operator between the two processes
        // (normalize() additionally inserts an ORIGIN for the dangling `in` channel)
        vFoo.type == 'PROCESS'
        vMap.type == 'OPERATOR'
        vBar.type == 'PROCESS'
        and:
        // the operator vertex is not dropped/collapsed
        result.vertices.count { it.type == 'OPERATOR' } == 1
        and:
        // ids are stable strings and unique, edges join via those ids
        result.vertices*.id.every { it instanceof String }
        result.vertices*.id.unique().size() == result.vertices.size()

        and:
        // the process -> operator -> process chain is intact
        result.edges.find { it.source == vFoo.id && it.target == vMap.id }
        result.edges.find { it.source == vMap.id && it.target == vBar.id }
    }

    def 'should expose fully-qualified process name and enclosing scope' () {
        given:
        def dag = new DAG()
        and:
        // a process invoked inside a subworkflow gets a scope-prefixed FQ name as its label,
        // while `workflow` holds the enclosing (sub)workflow scope
        def v = dag.createVertex(DAG.Type.PROCESS, 'RNASEQ:ALIGN:STAR')
        v.workflow = 'RNASEQ:ALIGN'

        when:
        def result = DagSerializer.toMap(dag)

        then:
        result.vertices.size() == 1
        result.vertices[0].label == 'RNASEQ:ALIGN:STAR'
        result.vertices[0].scope == 'RNASEQ:ALIGN'
    }

    def 'should fill in boundary vertices for dangling edges' () {
        given:
        def ch1 = Mock(DataflowChannel)
        def ch2 = Mock(DataflowChannel)
        and:
        def dag = new DAG()
        // a single process with one dangling input and one dangling output:
        // normalize() should create an ORIGIN source and a NODE termination
        dag.addVertex(
                DAG.Type.PROCESS,
                'foo',
                [ new DAG.ChannelHandler(channel: ch1) ],
                [ new DAG.ChannelHandler(channel: ch2) ] )

        when:
        def result = DagSerializer.toMap(dag)

        then:
        result.vertices*.type.toSet() == ['ORIGIN', 'PROCESS', 'NODE'] as Set
        and:
        // both edges now have a resolved source and target
        result.edges.every { it.source != null && it.target != null }
    }

}
