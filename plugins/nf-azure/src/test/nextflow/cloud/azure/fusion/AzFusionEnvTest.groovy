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

package nextflow.cloud.azure.fusion

import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.fusion.FusionConfig
import spock.lang.Specification

/**
 *
 * @author Alberto Miranda <alberto.miranda@seqera.io>
 */
class AzFusionEnvTest extends Specification {

    def setup() {
        SysEnv.push([:])  // <-- clear the system host env
    }

    def cleanup() {
        SysEnv.pop()      // <-- restore the system host env
    }

    def 'should return empty env'() {
        given:
        def provider = new AzFusionEnv()
        when:
        def env = provider.getEnvironment('aws', Mock(FusionConfig))
        then:
        env == Collections.emptyMap()
    }

    def 'should return env environment'() {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [azure: [storage: [accountName: 'x1']]]
        }
        and:

        when:
        def config = Mock(FusionConfig)
        def env = new AzFusionEnv().getEnvironment('az', config)
        then:
        env == [AZURE_STORAGE_ACCOUNT: 'x1']

    }

    def 'should return env environment with SAS token config'() {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [azure: [storage: [accountName: 'x1', sasToken: 'y1']]]
        }
        and:

        when:
        def config = Mock(FusionConfig)
        def env = new AzFusionEnv().getEnvironment('az', config)
        then:
        env == [AZURE_STORAGE_ACCOUNT: 'x1', AZURE_STORAGE_SAS_TOKEN: 'y1']

    }

    def 'should throw an exception when missing Azure Storage account name'() {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [azure: [storage: [sasToken: 'y1']]]
        }
        when:
        def config = Mock(FusionConfig)
        def env = new AzFusionEnv().getEnvironment('az', Mock(FusionConfig))
        then:
        thrown(IllegalArgumentException)
    }

    def 'should throw an exception when both account key and SAS token are present'() {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [azure: [storage: [accountName: 'x1', accountKey: 'y1', sasToken: 'z1']]]
        }
        when:
        def config = Mock(FusionConfig)
        def env = new AzFusionEnv().getEnvironment('az', Mock(FusionConfig))
        then:
        thrown(IllegalArgumentException)
    }
}
