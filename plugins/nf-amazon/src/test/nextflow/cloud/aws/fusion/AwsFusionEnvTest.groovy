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

package nextflow.cloud.aws.fusion

import nextflow.Global
import nextflow.SysEnv
import nextflow.fusion.FusionConfig
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsFusionEnvTest extends Specification {

    def setup() {
        Global.config = Collections.emptyMap()
    }

    def 'should return empty env' () {
        given:
        def provider = new AwsFusionEnv()
        when:
        def env = provider.getEnvironment('az', Mock(FusionConfig))
        then:
        env == Collections.emptyMap()
    }

    def 'should return env environment' () {
        given:
        SysEnv.push([AWS_ACCESS_KEY_ID: 'x1', AWS_SECRET_ACCESS_KEY: 'y1', AWS_S3_ENDPOINT: 'http://my-host.com'])
        and:

        when:
        def config = Mock(FusionConfig)
        def env = new AwsFusionEnv().getEnvironment('s3', Mock(FusionConfig))
        then:
        env == [AWS_S3_ENDPOINT:'http://my-host.com']

        when:
        config = Mock(FusionConfig) { exportStorageCredentials() >> true }
        env = new AwsFusionEnv().getEnvironment('s3', config)
        then:
        env == [AWS_ACCESS_KEY_ID: 'x1',
                AWS_SECRET_ACCESS_KEY: 'y1',
                AWS_S3_ENDPOINT:'http://my-host.com']

        cleanup:
        SysEnv.pop()
    }

    def 'should return env environment with SSE config' () {
        given:
        Global.config = [aws:[client: [storageEncryption:'aws:kms', storageKmsKeyId: 'xyz']]]
        and:

        when:
        def config = Mock(FusionConfig)
        def env = new AwsFusionEnv().getEnvironment('s3', Mock(FusionConfig))
        then:
        env == [FUSION_AWS_SERVER_SIDE_ENCRYPTION:'aws:kms', FUSION_AWS_SSEKMS_KEY_ID:'xyz']

        cleanup:
        Global.config = null
    }

    def 'should return env environment with session token' () {
        given:
        SysEnv.push([AWS_ACCESS_KEY_ID: 'x1', AWS_SECRET_ACCESS_KEY: 'y1', AWS_S3_ENDPOINT: 'http://my-host.com', AWS_SESSION_TOKEN: 'z1'])
        and:

        when:
        def config = Mock(FusionConfig)
        def env = new AwsFusionEnv().getEnvironment('s3', Mock(FusionConfig))
        then:
        env == [AWS_S3_ENDPOINT:'http://my-host.com']

        when:
        config = Mock(FusionConfig) { exportStorageCredentials() >> true }
        env = new AwsFusionEnv().getEnvironment('s3', config)
        then:
        env == [AWS_ACCESS_KEY_ID: 'x1',
                AWS_SECRET_ACCESS_KEY: 'y1',
                AWS_S3_ENDPOINT:'http://my-host.com',
                AWS_SESSION_TOKEN: 'z1']

        cleanup:
        SysEnv.pop()
    }
}
