/*
 * Copyright 2020-2021, Seqera Labs
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
 *
 */

package nextflow.extension

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.Session
import nextflow.plugin.Scoped
import spock.lang.Specification
import test.MockScriptRunner
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ChannelFactoryImplTest extends Specification {

    @Scoped('foo')
    class Ext1 extends ChannelExtensionPoint {

        int initCount
        Session initSession

        @Override
        void init(Session sess) {
            initSession = sess
            initCount += 1
        }

        DataflowWriteChannel alpha(List args) {
            def ch = new DataflowQueue()
            return CH.emitAndClose(ch, args)
        }

        DataflowWriteChannel omega() {
            Channel.value('Hello world')
        }

        DataflowWriteChannel plusOne( DataflowReadChannel ch ) {
            ch.map { it +1 }
        }
    }

    @Scoped('bar')
    class Ext2 extends ChannelExtensionPoint {
        int initCount
        Session initSession

        @Override
        void init(Session sess) {
            initSession = sess
            initCount += 1
        }

        DataflowWriteChannel omega(List args) {
            def ch = new DataflowQueue()
            return CH.emitAndClose(ch, args)
        }
    }

    def 'should invoke custom plugin factory' () {
        given:
        OperatorEx.channelFactories.clear()
        and:
        def ext1 = new Ext1(); def ext2 = new Ext2()
        ChannelFactoryImpl.allExtensions = [ext1,ext2]
        and:
        def SCRIPT = '''
        Channel.foo.alpha(['one','two','three'])
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        result.val == 'one'
        result.val == 'two'
        result.val == 'three'
        result.val == Channel.STOP
        and:
        ext1.initCount == 1
        ext1.initSession instanceof Session
        and:
        ext2.initCount == 0
        ext2.initSession == null
        cleanup:
        ChannelFactoryImpl.allExtensions=null
    }


    def 'should invoke multiple extensions' () {
        given:
        OperatorEx.channelFactories.clear()
        and:
        def ext1 = new Ext1(); def ext2 = new Ext2()
        ChannelFactoryImpl.allExtensions = [ext1,ext2]
        and:
        def SCRIPT = '''
            def ch1 = channel.foo.alpha([1,2,3])
            def ch2 = channel.bar.omega(['X','Y','Z'])
            
            process sayHello {
              input:
                val x from ch1
                val y from ch2
              output: 
                val z into ch3 
              exec:
                z = "$x $y"
            }
            
            ch3.toSortedList()
            '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        result.val == ['1 X', '2 Y', '3 Z']

        and:
        ext1.initCount == 1
        ext1.initSession instanceof Session
        and:
        ext2.initCount == 1
        ext2.initSession instanceof Session
        cleanup:
        ChannelFactoryImpl.allExtensions=null

    }

    def 'should invoke operator extension' () {
        given:
        OperatorEx.channelFactories.clear()
        and:
        def ext1 = new Ext1();
        ChannelFactoryImpl.allExtensions = [ext1]
        and:
        def SCRIPT = '''
            channel
                .fromList([1,2,3])
                .plusOne()
            '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        result.val == 2
        result.val == 3
        result.val == 4

        and:
        ext1.initCount == 1
        ext1.initSession instanceof Session

        cleanup:
        ChannelFactoryImpl.allExtensions=null
    }

}
