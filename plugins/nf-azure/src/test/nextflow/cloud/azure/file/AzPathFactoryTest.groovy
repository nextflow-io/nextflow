/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cloud.azure.file

import nextflow.Global
import nextflow.Session
import nextflow.cloud.azure.nio.AzPath
import spock.lang.Requires
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Requires({System.getenv('AZURE_STORAGE_ACCOUNT_NAME') && System.getenv('AZURE_STORAGE_ACCOUNT_KEY')})
class AzPathFactoryTest extends Specification {

    def 'should create az azure path' () {
        given:
        def CONFIG = [azure: [
                storage: [
                        accountKey: System.getenv('AZURE_STORAGE_ACCOUNT_KEY'),
                        accountName: System.getenv('AZURE_STORAGE_ACCOUNT_NAME'),
                ]
        ]]
        Global.session = Mock(Session) { getConfig() >> CONFIG }
        and:

        when:
        def path = AzPathFactory.parse(AZ_URI)
        then:
        path instanceof AzPath
        (path as AzPath).containerName == CONTAINER
        (path as AzPath).blobName() == BLOB

        when:
        def ret = AzPathFactory.getUriString(path)
        then:
        ret == AZ_URI

        cleanup:
        Global.session = null

        where:
        AZ_URI                          | CONTAINER     | BLOB
        'az://my-data/foo/bar'          | 'my-data'     | 'foo/bar'
        'az://my-data/data/*{1,2}.fq.gz'| 'my-data'     | 'data/*{1,2}.fq.gz'
    }



    def 'should throw illegal path' () {
        given:
        def CONFIG = [azure: [
                storage: [
                        accountKey: System.getenv('AZURE_STORAGE_ACCOUNT_KEY'),
                        accountName: System.getenv('AZURE_STORAGE_ACCOUNT_NAME'),
                ]
        ]]
        Global.session = Mock(Session) { getConfig() >> CONFIG }
        and:

        when:
        AzPathFactory.parse('az:///fooo')
        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Invalid Azure path URI - make sure the schema prefix does not container more than two slash characters - offending value: az:///fooo'
    }
}
