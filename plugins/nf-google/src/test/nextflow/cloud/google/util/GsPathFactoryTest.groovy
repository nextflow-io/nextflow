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
 */

package nextflow.cloud.google.util

import com.google.cloud.storage.StorageOptions
import nextflow.Global
import nextflow.Session
import nextflow.cloud.google.GoogleOpts
import nextflow.cloud.google.config.GoogleRetryOpts
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GsPathFactoryTest extends Specification {

    @Unroll
    def 'should create gs path #PATH' () {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [google:[project:'foo', region:'x']]
        }
        and:
        def factory = new GsPathFactory()

        expect:
        factory.parseUri(PATH).toUriString() == PATH
        factory.parseUri(PATH).toString() == STR

        where:
        _ | PATH                | STR
        _ | 'gs://foo'          | ''
        _ | 'gs://foo/bar'      | '/bar'
        _ | 'gs://foo/bar/'     | '/bar/'   // <-- bug or feature ?
        _ | 'gs://foo/b a r'    | '/b a r'
        _ | 'gs://f o o/bar'    | '/bar'
        _ | 'gs://f_o_o/bar'    | '/bar'
    }

    def 'should use requester pays' () {
        given:
        def session = Mock(Session) {
            getConfig() >> [google:[project:'foo', region:'x', enableRequesterPaysBuckets:true]]
        }
        and:
        def opts = GoogleOpts.fromSession(session)

        when:
        def storageConfig = GsPathFactory.getCloudStorageConfig(opts)

        then:
        storageConfig.userProject() == 'foo'
    }

    def 'should not use requester pays' () {
        given:
        def session = new Session()
        session.config = [google:[project:'foo', region:'x', lifeSciences: [:]]]
        and:
        def opts = GoogleOpts.fromSession(session)

        when:
        def storageConfig = GsPathFactory.getCloudStorageConfig(opts)

        then:
        storageConfig.userProject() == null
    }

    def 'should apply http timeout settings from config' () {
        given:
        def session = Mock(Session) {
            getConfig() >> [google:[httpConnectTimeout: CONNECT, httpReadTimeout: READ]]
        }
        and:
        def policy = new GoogleRetryOpts([:])
        def retrySettings = StorageOptions.getDefaultRetrySettings()
            .toBuilder()
            .setMaxAttempts(policy.maxAttempts)
            .setRetryDelayMultiplier(policy.multiplier)
            .setTotalTimeout(org.threeten.bp.Duration.ofSeconds(policy.maxDelaySecs()))
            .build()
        and:
        def opts = GoogleOpts.fromSession(session)
        and:
        def storageOptions = GsPathFactory.getCloudStorageOptions(opts)
        and:
        def transportOptions = StorageOptions.getDefaultHttpTransportOptions().toBuilder()
        if( CONNECT ) transportOptions.setConnectTimeout( CONNECT_MILLIS )
        if( READ ) transportOptions.setReadTimeout( READ_MILLIS )

        expect:
        storageOptions == StorageOptions.getDefaultInstance().toBuilder()
            .setTransportOptions(transportOptions.build())
            .setRetrySettings(retrySettings)
            .build()

        where:
        CONNECT | CONNECT_MILLIS | READ  | READ_MILLIS
        null    | 60000          | null  | 60000
        '30s'   | 30000          | '30s' | 30000
        '60s'   | 60000          | '60s' | 60000
    }

    def 'should apply retry settings' () {
        given:
        def session = Mock(Session) {
            getConfig() >> [google:[storage:[retryPolicy: [maxAttempts: 5, maxDelay:'50s', multiplier: 500]]]]
        }

        when:
        def opts = GoogleOpts.fromSession(session)
        then:
        opts.storageOpts.retryPolicy.maxAttempts == 5
        opts.storageOpts.retryPolicy.maxDelaySecs() == 50
        opts.storageOpts.retryPolicy.multiplier == 500d

    }
}
