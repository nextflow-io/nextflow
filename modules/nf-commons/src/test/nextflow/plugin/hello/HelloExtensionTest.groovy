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

import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.plugin.hello.HelloExtension
import spock.lang.Specification


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
class HelloExtensionTest extends Specification {

    def "should create a channel from hello"(){
        given:
        def session = Mock(Session)
        Global.session = session

        and:
        def helloExtension = new HelloExtension(); helloExtension.init(session)

        when:
        def result = helloExtension.reverse("Hi")

        then:
        result.val == 'iH'
        result.val == Channel.STOP
    }

    def "should consume a message from script"(){

        given:
        def session = Mock(Session)
        Global.session = session

        and:
        def helloExtension = new HelloExtension(); helloExtension.init(session)

        and:
        def ch = new DataflowQueue()
        ch.bind('Goodbye folks')
        ch.bind( Channel.STOP )

        when:
        def result = helloExtension.goodbye(ch)

        then:
        result.val == 'Goodbye folks'
        result.val == Channel.STOP
    }
}
