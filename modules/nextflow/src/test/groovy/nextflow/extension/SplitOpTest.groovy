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

package nextflow.extension

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.splitter.FastqSplitter
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SplitOpTest extends Specification {

    static abstract class FakeDataflowQueue implements DataflowReadChannel, DataflowWriteChannel { }

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
        def METHOD = 'splitFastq'
        def SOURCE = Mock(DataflowReadChannel)
        def params = [:]
        def op = Spy(SplitOp, constructorArgs:[SOURCE, METHOD, params])

        def splitter = Mock(FastqSplitter)
        def OUTPUT = Mock(DataflowWriteChannel)

        when:
        op.apply()
        then:
        1 * op.getOrCreateDataflowQueue([autoClose:false]) >> OUTPUT
        1 * op.createSplitter(METHOD, [autoClose: false, into:OUTPUT]) >> splitter
        1 * op.applySplittingOperator(SOURCE, OUTPUT, splitter)
        0 * splitter.setMultiSplit(_) >> null
        0 * splitter.setEmitSplitIndex(_) >> null

    }

    def 'should invoke paired-end split' () {
        given:
        def METHOD = 'splitFastq'
        def SOURCE = Mock(DataflowReadChannel)
        def params = [pe:true]
        def op = Spy(SplitOp, constructorArgs:[SOURCE, METHOD, params])

        def copy1 = Mock(FakeDataflowQueue)
        def copy2 = Mock(FakeDataflowQueue)
        def out1 = Mock(DataflowWriteChannel)
        def out2 = Mock(DataflowWriteChannel)
        def splitter1 = Mock(FastqSplitter)
        def splitter2 = Mock(FastqSplitter)

        when:
        def result = op.apply()
        then:
        1 * op.splitMultiEntries()
        1 * op.createSourceCopies(SOURCE,2) >> [copy1, copy2]

        1 * op.splitSingleEntry(copy1, [elem:-1, autoClose: false])
        1 * op.splitSingleEntry(copy2, [elem:-2, autoClose: false])

        1 * op.getOrCreateDataflowQueue( _ as Map ) >> out1
        1 * op.getOrCreateDataflowQueue( _ as Map ) >> out2

        1 * op.createSplitter(METHOD, [elem:-1, autoClose: false, into:out1] ) >> splitter1
        1 * op.createSplitter(METHOD, [elem:-2, autoClose: false, into:out2] ) >> splitter2

        1 * op.applySplittingOperator(copy1, out1, splitter1) >> null
        1 * op.applySplittingOperator(copy2, out2, splitter2) >> null

        1 * splitter1.setEmitSplitIndex(true)
        1 * splitter2.setEmitSplitIndex(true)

        1 * splitter1.setMultiSplit(true)
        1 * splitter2.setMultiSplit(true)
        1 * op.applyMergingOperator([out1, out2], _ as DataflowQueue, [-1,-2]) >> null
        then:
        result instanceof DataflowQueue
    }


}
