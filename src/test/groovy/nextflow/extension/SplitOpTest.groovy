package nextflow.extension

import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.splitter.AbstractSplitter
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SplitOpTest extends Specification {

    def 'should validate params' () {

        given:
        SplitOp op

        when:
        op = new SplitOp(Mock(DataflowReadChannel), 'splitFasta', [:])
        then:
        op.methodName == 'splitFasta'
        op.params == [autoClose:false]
        !op.pairedEnd
        !op.multiSplit
        op.indexes == null


        when:
        op = new SplitOp(Mock(DataflowReadChannel), 'splitFasta', [elem:1])
        then:
        op.methodName == 'splitFasta'
        !op.pairedEnd
        !op.multiSplit
        !op.indexes

        when:
        op = new SplitOp(Mock(DataflowReadChannel), 'splitFasta', [elem:[1,2,4]])
        then:
        op.methodName == 'splitFasta'
        !op.pairedEnd
        op.multiSplit
        op.indexes == [1,2,4]

        when:
        op = new SplitOp(Mock(DataflowReadChannel), 'splitFastq', [pe:true])
        then:
        op.methodName == 'splitFastq'
        op.pairedEnd
        op.multiSplit
        op.indexes == [-1,-2]

        when:
        new SplitOp(Mock(DataflowReadChannel), 'splitFasta', [autoClose: true])
        then:
        thrown(IllegalArgumentException)

        when:
        new SplitOp(Mock(DataflowReadChannel), 'splitFasta', [into: 'any'])
        then:
        thrown(IllegalArgumentException)

        when:
        new SplitOp(Mock(DataflowReadChannel), 'splitFastq', [pe:true, elem:1])
        then:
        thrown(IllegalArgumentException)

        when:
        new SplitOp(Mock(DataflowReadChannel), 'splitFasta', [pe:true])
        then:
        thrown(IllegalArgumentException)
    }


    def 'should invoke single split' () {

        given:
        def METHOD = 'splitFasta'
        def SOURCE = Mock(DataflowReadChannel)
        def params = [:]
        def op = Spy(SplitOp, constructorArgs:[SOURCE, METHOD, params])

        def SPLITTER = Mock(AbstractSplitter)
        def OUTPUT = Mock(DataflowWriteChannel)

        when:
        op.apply()
        then:
        1 * op.getOrCreateDataflowQueue([autoClose:false]) >> OUTPUT
        1 * op.createSplitter(METHOD, [autoClose: false]) >> SPLITTER
        1 * op.applySplittingOperator(SOURCE, OUTPUT, SPLITTER)


    }
}
