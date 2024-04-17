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

import com.amazonaws.AmazonClientException
import com.amazonaws.auth.AWSCredentials
import com.amazonaws.auth.AWSCredentialsProvider
import com.amazonaws.auth.AnonymousAWSCredentials
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 * AWS credentials provider that delegates the credentials to the
 * specified provider class and fallback to the {@link AnonymousAWSCredentials}
 * when no credentials are available.
 *
 * See also {@link com.amazonaws.services.s3.S3CredentialsProviderChain}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class S3CredentialsProvider implements AWSCredentialsProvider {

    private AWSCredentialsProvider target

    private volatile AWSCredentials anonymous

    S3CredentialsProvider(AWSCredentialsProvider target) {
        this.target = target
    }

    @Override
    AWSCredentials getCredentials() {
        if( anonymous!=null )
            return anonymous
        try {
            return target.getCredentials();
        } catch (AmazonClientException e) {
            log.debug("No AWS credentials available - falling back to anonymous access");
        }
        return anonymous=new AnonymousAWSCredentials()
    }

    @Override
    void refresh() {
        target.refresh()
    }
}
