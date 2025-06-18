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

import software.amazon.awssdk.auth.credentials.AwsCredentials
import software.amazon.awssdk.auth.credentials.AwsCredentialsProvider
import software.amazon.awssdk.auth.credentials.AnonymousCredentialsProvider
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 * AWS credentials provider that delegates the credentials to the
 * specified provider class and fallback to the {@link AnonymousCredentialsProvider}
 * when no credentials are available.
 *
 * See also {@link software.amazon.awssdk.auth.credentials.AwsCredentialsProviderChain}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class S3CredentialsProvider implements AwsCredentialsProvider {

    private AwsCredentialsProvider target

    private volatile AwsCredentials anonymous

    S3CredentialsProvider(AwsCredentialsProvider target) {
        this.target = target
    }

    @Override
    AwsCredentials resolveCredentials() {
        if (anonymous != null) {
            return anonymous
        }
        try {
            return target.resolveCredentials()
        } catch (Exception e) {
            log.debug("No AWS credentials available - falling back to anonymous access")
        }
        anonymous = AnonymousCredentialsProvider.create().resolveCredentials()
        return anonymous
    }

}
