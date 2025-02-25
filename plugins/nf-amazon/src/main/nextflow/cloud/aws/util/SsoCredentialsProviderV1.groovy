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

package nextflow.cloud.aws.util

import com.amazonaws.auth.AWSCredentials
import com.amazonaws.auth.AWSCredentialsProvider
import com.amazonaws.auth.BasicAWSCredentials
import com.amazonaws.auth.BasicSessionCredentials
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import software.amazon.awssdk.auth.credentials.AwsSessionCredentials
import software.amazon.awssdk.auth.credentials.ProfileCredentialsProvider

/**
 * Adapter for the SSO credentials provider from the SDK v2.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class SsoCredentialsProviderV1 implements AWSCredentialsProvider {

    private ProfileCredentialsProvider delegate

    SsoCredentialsProviderV1() {
        this.delegate = ProfileCredentialsProvider.create()
    }

    SsoCredentialsProviderV1(String profile) {
        this.delegate = ProfileCredentialsProvider.create(profile)
    }

    @Override
    AWSCredentials getCredentials() {
        final credentials = delegate.resolveCredentials()

        if( credentials instanceof AwsSessionCredentials )
            new BasicSessionCredentials(
                    credentials.accessKeyId(),
                    credentials.secretAccessKey(),
                    credentials.sessionToken())

        else
            new BasicAWSCredentials(
                    credentials.accessKeyId(),
                    credentials.secretAccessKey())
    }

    @Override
    void refresh() {
        throw new UnsupportedOperationException()
    }
}
