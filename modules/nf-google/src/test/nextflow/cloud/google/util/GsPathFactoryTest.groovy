/*
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

package nextflow.cloud.google.util

import nextflow.Global
import nextflow.Session
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GsPathFactoryTest extends Specification {

    def 'should create gs path' () {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [google:[project:'foo', region:'x']]
        }
        and:
        def factory = new GsPathFactory()

        expect:
        factory.parseUri(PATH).toUriString() == PATH

        where:
        _ | PATH
        _ | 'gs://foo'
        _ | 'gs://foo/bar'
        _ | 'gs://foo/b a r'
        _ | 'gs://f o o/bar'
        _ | 'gs://f_o_o/bar'
    }

    def 'should use requester pays' () {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [google:[project:'foo', region:'x', enableRequesterPaysBuckets:true]]
        }

        when:
        def storageConfig = GsPathFactory.getCloudStorageConfig()

        then:
        storageConfig.userProject() == 'foo'
    }

    def 'should not use requester pays' () {
        given:
        def sess = new Session()
        sess.config = [google:[project:'foo', region:'x', lifeSciences: [:]]]
        Global.session = sess

        when:
        def storageConfig = GsPathFactory.getCloudStorageConfig()

        then:
        storageConfig.userProject() == null
    }
}
