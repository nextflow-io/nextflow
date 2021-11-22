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

package nextflow.sraql

import nextflow.Channel
import nextflow.extension.ChannelExtensionDelegate
import spock.lang.Timeout
import test.BaseSpec
import test.MockScriptRunner

/**
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
@Timeout(20)
class SraqlDslTest extends BaseSpec {

    def setup() {
        ChannelExtensionDelegate.reloadExtensionPoints()
    }

    def 'should perform a query and create a channel' () {
        def config = [sraql: [source: 'google-bigquery']]

        when:
        def SCRIPT = '''
            def bioProjectId = 'PRJNA494931'
            def sraql = "SELECT *  FROM `nih-sra-datastore.sra.metadata` WHERE  bioproject='$bioProjectId';"
            channel.sraql.fromQuery(sraql)
            '''
        and:
        def result = new MockScriptRunner(config).setScript(SCRIPT).execute()

        then:
        result.val.acc == 'SRR7974377'
        result.val.acc == 'SRR7974375'
        result.val.acc == 'SRR7974376'
        result.val == Channel.STOP
    }

}
