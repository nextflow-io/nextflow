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

package nextflow.plugin.hello

import nextflow.Channel
import nextflow.plugin.extension.PluginExtensionProvider
import spock.lang.Specification
import spock.lang.Timeout
import test.MockScriptRunner


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
@Timeout(10)
class HelloDslTest extends Specification{

    def setup () {
        def ext = new PluginExtensionProvider(){
            def loadForTest(){
                install()
                loadPluginExtensionMethods("nf-foo", new HelloExtension(), [reverse:'reverse', goodbye:'goodbye'])
            }
        }
        ext.loadForTest()
    }

    def 'should perform a hi and create a channel' () {
        when:
        def SCRIPT = '''
            channel.reverse('hi!') 
            '''
        and:
        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
        then:
        result.val == '!ih'
        result.val == Channel.STOP
    }

    def 'should store a goodbye' () {
        when:
        def SCRIPT = '''
            channel
                .of('Bye bye folks')
                .goodbye() 
            '''
        and:
        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
        then:
        result.val == 'Bye bye folks'
        result.val == Channel.STOP
        
    }
}
