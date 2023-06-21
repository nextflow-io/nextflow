/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.cloud.aws.util

import com.amazonaws.auth.AWSCredentials
import com.amazonaws.auth.AWSCredentialsProvider
import com.amazonaws.auth.BasicAWSCredentials
import com.amazonaws.auth.BasicSessionCredentials
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import software.amazon.awssdk.auth.credentials.AwsSessionCredentials
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.sso.SsoClient
import software.amazon.awssdk.services.sso.auth.SsoCredentialsProvider

/**
 * Adapter for the SSO credentials provider from the SDK v2.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class SsoCredentialsProviderV1 implements AWSCredentialsProvider {

    private SsoCredentialsProvider delegate

    SsoCredentialsProviderV1(String region) {
        final ssoClient = SsoClient.builder()
                .region( Region.of(region) )
                .build()

        this.delegate = SsoCredentialsProvider.builder()
                .ssoClient( ssoClient )
                .build()
    }

    @Override
    AWSCredentials getCredentials() {
        final credentials = delegate.resolveCredentials()

        if( credentials instanceof AwsSessionCredentials ) {
            final sessionCredentials = (AwsSessionCredentials) credentials
            new BasicSessionCredentials(
                    sessionCredentials.accessKeyId(),
                    sessionCredentials.secretAccessKey(),
                    sessionCredentials.sessionToken())
        }
        else
            new BasicAWSCredentials(
                    credentials.accessKeyId(),
                    credentials.secretAccessKey())
    }

    @Override
    void refresh() {}
}
