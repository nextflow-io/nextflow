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

import com.google.cloud.bigquery.BigQueryOptions
import com.google.cloud.bigquery.JobId
import com.google.cloud.bigquery.JobInfo
import com.google.cloud.bigquery.QueryJobConfiguration
import groovyx.gpars.dataflow.DataflowQueue
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
        def bigquery = BigQueryOptions.getDefaultInstance().getService()
        def queryConfig = QueryJobConfiguration.newBuilder(
                "SELECT * "
                        + " FROM `nih-sra-datastore.sra.metadata` as s "
                        + " WHERE s.organism = 'Mycobacterium tuberculosis'"
                        + " AND s.consent='public' "
                        + " AND s.sra_study='ERP124850' "
                        + " LIMIT 3"
        )
                .setUseLegacySql(false)
                .build()

        def jobId = JobId.of(UUID.randomUUID().toString());
        def queryJob = bigquery
                .create(JobInfo.newBuilder(queryConfig).setJobId(jobId).build())
                .waitFor()

        and:
        def session = Mock(Session) {
            getConfig() >> [sraql: [source: 'google-bigquery']]
        }
        def sraqlExtension = new ChannelSraqlExtension(); sraqlExtension.init(session)

        when:
        def result = new DataflowQueue()
        def query = "SELECT *  FROM `nih-sra-datastore.sra.metadata` WHERE organism = 'Mycobacterium tuberculosis' LIMIT 3;"
        new QueryHandler()
                .withTarget(result)
                .withStatement(query)
                .perform()

        then:
        result.val[0] == ['acc', 'ERR4796597']
        result.val[0] == ['acc', 'ERR4797168']
        result.val[0] == ['acc', 'ERR4797173']
        result.val == Channel.STOP

    }
}
