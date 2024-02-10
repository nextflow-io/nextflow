/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.extension.plugin

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.Global
import nextflow.NextflowMeta
import nextflow.Session
import nextflow.extension.CH
import nextflow.plugin.extension.PluginExtensionPoint
import nextflow.plugin.extension.PluginExtensionProvider
import spock.lang.Specification
import test.MockScriptRunner
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ChannelFactoryInstanceTest extends Specification {

    class Ext1 extends PluginExtensionPoint {

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
            def result = new DataflowQueue()
            result.bind(2)
            result.bind(3)
            result.bind(4)
            result.bind(Channel.STOP)
            return result
        }
    }

    class Ext2 extends PluginExtensionPoint {
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
        Global.session = Mock(Session)
        and:
        def ext1 = new Ext1(); def ext2 = new Ext2()
        new PluginExtensionProvider()
                .install()
                .loadPluginExtensionMethods("nf-foo",ext1, ['alpha':'alpha'])
        and:
        def SCRIPT = '''
        Channel.alpha(['one','two','three'])
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
        PluginExtensionProvider.reset()
    }

    def 'should invoke alias in custom plugin factory' () {
        given:
        def ext1 = new Ext1(); def ext2 = new Ext2()
        new PluginExtensionProvider()
                .install()
                .loadPluginExtensionMethods("nf-foo",ext1, ['alpha':'thisIsAnAliasToAlpha'])
        and:
        def SCRIPT = '''
        Channel.thisIsAnAliasToAlpha(['one','two','three'])
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
        PluginExtensionProvider.reset()
    }

    def 'should invoke multiple extensions' () {
        given:
        NextflowMeta.instance.enableDsl2()
        Global.session = Mock(Session)
        and:
        def ext1 = new Ext1(); def ext2 = new Ext2()
        new PluginExtensionProvider()
                .install()
                .loadPluginExtensionMethods("nf-foo",ext1, ['alpha':'alpha'])
                .loadPluginExtensionMethods("nf-foo",ext2, ['omega':'omega'])
        and:
        def SCRIPT = '''
            def ch1 = channel.alpha([1,2,3])
            def ch2 = channel.omega(['X','Y','Z'])
            
            process sayHello {
              input:
                val x
                val y
              output: 
                val z
              exec:
                z = "$x $y"
            }
            
            workflow {
              main: sayHello(ch1, ch2)
              emit: sayHello.out.toSortedList()
            }
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
        PluginExtensionProvider.reset()
        NextflowMeta.instance.disableDsl2()
    }

    def 'should invoke operator extension' () {
        given:
        def ext1 = new Ext1();
        new PluginExtensionProvider()
                .install()
                .loadPluginExtensionMethods("nf-foo",ext1, ['plusOne':'plusOne'])
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
        PluginExtensionProvider.reset()
    }

}
