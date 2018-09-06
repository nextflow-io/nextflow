/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
