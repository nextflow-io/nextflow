package nextflow.extension

import spock.lang.Specification

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class OpCallTest extends Specification {

    static class MyExt extends OperatorEx{

        DataflowWriteChannel foo1(DataflowReadChannel source) {
            new DataflowQueue() << 1 << 2
        }

        DataflowWriteChannel foo2(DataflowReadChannel source, DataflowReadChannel channel2) {
            new DataflowQueue() << 1 << 2
        }

        String foo3(DataflowReadChannel source, List<DataflowReadChannel> channel3) {
            return 'hello'
        }

        String foo4(DataflowReadChannel source, DataflowReadChannel[] channel4) {
            return 'hello'
        }

        DataflowWriteChannel foo5(DataflowReadChannel source, DataflowWriteChannel channel2) {
            new DataflowQueue() << 1 << 2
        }

        List<DataflowReadChannel> foo6(DataflowReadChannel source, List<DataflowWriteChannel> channel3) {
            [ new DataflowQueue() ]
        }

        List<DataflowWriteChannel> foo7(DataflowReadChannel source, DataflowWriteChannel[] channel4) {
            [ new DataflowQueue() ]
        }
    }

    def 'should invoke and capture in/out definitions 1' () {

        given:
        def ext = new MyExt()
        def ch1 = new DataflowVariable()
        def op = new OpCall(ext, ch1, 'foo1', [] as Object[])

        when:
        def result = op.call()
        then:
        op.inputs == [ch1] as Set
        op.outputs == [result] as Set
        result instanceof DataflowQueue
        result.val == 1
        result.val == 2 

    }

    def 'should invoke and capture in/out definitions 2' () {

        given:
        def ext = new MyExt()
        def ch1 = new DataflowVariable()
        def ch2 = new DataflowVariable()
        def op = new OpCall(ext, ch1, 'foo2', [ch2] as Object[])

        when:
        def result = op.call()
        then:
        op.inputs == [ch1,ch2] as Set
        op.outputs == [result] as Set
        result instanceof DataflowQueue
        result.val == 1
        result.val == 2
    }

    def 'should invoke and capture in/out definitions 3' () {

        given:
        def ext = new MyExt()
        def ch1 = new DataflowVariable()
        def ch2 = new DataflowVariable()
        def ch3 = new DataflowVariable()
        def ch4 = new DataflowVariable()
        def op = new OpCall(ext, ch1, 'foo3', [[ch2,ch3,ch4]] as Object[])

        when:
        def result = op.call()
        then:
        op.inputs == [ch1,ch2,ch3,ch4] as Set
        op.outputs == [] as Set
        result == "hello"

    }

    def 'should invoke and capture in/out definitions 4' () {

        given:
        def ext = new MyExt()
        def ch1 = new DataflowVariable()
        def ch2 = new DataflowVariable()
        def ch3 = new DataflowVariable()
        def ch4 = new DataflowVariable()
        def op = new OpCall(ext, ch1, 'foo4', [[ch2,ch3,ch4] as DataflowReadChannel[]] as Object[])

        when:
        def result = op.call()
        then:
        op.inputs == [ch1,ch2,ch3,ch4] as Set
        op.outputs == [] as Set
        result == "hello"

    }

    def 'should invoke and capture in/out definitions 5' () {

        given:
        def ext = new MyExt()
        def ch1 = new DataflowVariable()
        def ch2 = new DataflowVariable()
        def op = new OpCall(ext, ch1, 'foo5', [ch2] as Object[])

        when:
        def result = op.call()
        then:
        op.inputs == [ch1] as Set
        op.outputs == [ch2, result] as Set

    }

    def 'should invoke and capture in/out definitions 6' () {

        given:
        def ext = new MyExt()
        def ch1 = new DataflowVariable()
        def ch2 = new DataflowVariable()
        def ch3 = new DataflowVariable()
        def op = new OpCall(ext, ch1, 'foo6', [ [ch2, ch3]] as Object[])

        when:
        def result = op.call()
        then:
        op.inputs == [ch1] as Set
        op.outputs == [ch2, ch3] as Set

    }

    def 'should invoke and capture in/out definitions 7' () {

        given:
        def ext = new MyExt()
        def ch1 = new DataflowVariable()
        def ch2 = new DataflowVariable()
        def ch3 = new DataflowVariable()
        def op = new OpCall(ext, ch1, 'foo7', [ [ch2, ch3] as DataflowWriteChannel[]] as Object[])

        when:
        def result = op.call()
        then:
        op.inputs == [ch1] as Set
        op.outputs == [ch2, ch3, result[0]] as Set

    }


    static class Dummy {

        void m1(DataflowReadChannel ch1 )  { }
        void m2(DataflowReadChannel ch1, DataflowReadChannel ch2)  { }
        void m3(DataflowReadChannel ch1, String foo ) { }
        void m4(DataflowReadChannel ch1, List<DataflowReadChannel> many ) { }
        void m5(DataflowReadChannel ch1, DataflowReadChannel... many ) { }
        void m6(DataflowReadChannel ch1, List<String> list ) { }
        void m7(DataflowReadChannel ch1, List list ) { }

        void r1() { }
        DataflowReadChannel r2() { }
        DataflowWriteChannel r3 () { }
        List<DataflowWriteChannel> r4() { }
        DataflowWriteChannel[] r5() { }
        List r6() { }


    }

    def 'should get param type' () {

        given:
        def ext = new OperatorEx()
        def op = new OpCall(ext, Mock(DataflowReadChannel), 'foo', Object[])

        expect:
        !op.declaresParamType(DataflowReadChannel, Dummy.getMethod('m1', [DataflowReadChannel] as Class[]))
        op.declaresParamType(DataflowReadChannel, Dummy.getMethod('m2', [DataflowReadChannel,DataflowReadChannel] as Class[]))
        !op.declaresParamType(DataflowReadChannel, Dummy.getMethod('m3', [DataflowReadChannel,String] as Class[]))
        //op.declaresParamType(DataflowReadChannel, Dummy.getMethod('m4', [DataflowReadChannel,LIST_DATAFLOW_READ.class] as Class[]))
        !op.declaresParamType(DataflowReadChannel, Dummy.getMethod('m6', [DataflowReadChannel,List] as Class[]))
        !op.declaresParamType(DataflowReadChannel, Dummy.getMethod('m7', [DataflowReadChannel,List] as Class[]))

    }

    def 'should get return type' () {

        given:
        def ext = new OperatorEx()
        def op = new OpCall(ext, Mock(DataflowReadChannel), 'foo', Object[])

        expect:
        !op.declaresReturnType(DataflowReadChannel, Dummy.getMethod('r1', [] as Class[]))
        op.declaresReturnType(DataflowReadChannel, Dummy.getMethod('r2', [] as Class[]))
        op.declaresReturnType(DataflowWriteChannel, Dummy.getMethod('r3', [] as Class[]))
        op.declaresReturnType(DataflowWriteChannel, Dummy.getMethod('r4', [] as Class[]))
        !op.declaresReturnType(DataflowReadChannel, Dummy.getMethod('r4', [] as Class[]))
        op.declaresReturnType(DataflowWriteChannel, Dummy.getMethod('r5', [] as Class[]))
        !op.declaresReturnType(DataflowWriteChannel, Dummy.getMethod('r6', [] as Class[]))
    }
}
