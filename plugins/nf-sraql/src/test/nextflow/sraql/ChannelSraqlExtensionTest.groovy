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
import nextflow.Session
import spock.lang.Specification
import spock.lang.Timeout
/**
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
@Timeout(20)
class ChannelSraqlExtensionTest extends Specification {

    def 'should read the config for data source, execute query and create channel from query'() {
        given:
        def session = Mock(Session) {
            getConfig() >> [sraql: [source: 'google-bigquery']]
        }
        def sraqlExtension = new ChannelSraqlExtension(); sraqlExtension.init(session)


        and:
        def queryString = "SELECT *  FROM `nih-sra-datastore.sra.metadata` WHERE  bioproject='PRJNA494931';"

        when:
        def result = sraqlExtension.fromQuery(queryString)

        then:
        result.val.acc == 'SRR7974377'
        result.val.acc == 'SRR7974375'
        result.val.acc == 'SRR7974376'
        result.val == Channel.STOP
    }
}
