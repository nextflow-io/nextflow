/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.cloud.aws.batch

import com.amazonaws.services.batch.AWSBatchClient
import com.amazonaws.services.batch.model.DescribeJobDefinitionsRequest
import com.amazonaws.services.batch.model.DescribeJobDefinitionsResult
import com.amazonaws.services.batch.model.DescribeJobsRequest
import com.amazonaws.services.batch.model.DescribeJobsResult
import nextflow.util.ThrottlingExecutor
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsBatchProxyTest extends Specification {

    def 'should get client instance' () {

        given:
        def client = Mock(AWSBatchClient)
        def exec = Mock(ThrottlingExecutor)

        when:
        def c = new AwsBatchProxy(client,exec).client
        then:
        0 * exec._
        c == client

        when:
        def d = new AwsBatchProxy(client,exec).getClient()
        then:
        0 * exec._
        d == client

    }

    def 'should invoke executor with normal priority' () {

        given:
        def client = Mock(AWSBatchClient)
        def exec = Mock(ThrottlingExecutor)
        def req = Mock(DescribeJobDefinitionsRequest)
        def resp = Mock(DescribeJobDefinitionsResult)
        def ZERO = 0 as byte

        when:
        def result = new AwsBatchProxy(client,exec).describeJobDefinitions(req)
        then:
        1 * exec.doInvoke1(client, 'describeJobDefinitions', [req] as Object[], ZERO) >> resp

        result == resp

    }

    def 'should invoke executor with higher priority' () {

        given:
        def client = Mock(AWSBatchClient)
        def exec = Mock(ThrottlingExecutor)
        def req = Mock(DescribeJobsRequest)
        def resp = Mock(DescribeJobsResult)
        def _10 = 10 as byte

        when:
        def result = new AwsBatchProxy(client,exec).describeJobs(req)
        then:
        1 * exec.doInvoke1(client, 'describeJobs', [req] as Object[], _10) >> resp

        result == resp

    }

}
