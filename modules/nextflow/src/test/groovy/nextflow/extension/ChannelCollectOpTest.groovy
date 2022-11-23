/*
 * Copyright 2020-2022, Seqera Labs
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

import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.exception.ScriptCompilationException
import org.junit.Rule
import spock.lang.Ignore
import test.Dsl2Spec
import test.OutputCapture

/**
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class ChannelCollectOpTest extends Dsl2Spec  {

    @Rule
    OutputCapture capture = new OutputCapture()


    /**
     * test chaining collect method
     */

    def 'test collect' () {

        when:
        def result = dsl_eval'''            
            Channel.of('a','b') | collect { it.toUpperCase() } | view
            '''
        then:
        result.val == ['A', 'B']
    }

}
