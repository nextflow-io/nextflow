/*
 * Copyright 2013-2025, Seqera Labs
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
package nextflow.cloud.aws.nio.util

import nextflow.cloud.aws.config.AwsConfig
import software.amazon.awssdk.auth.signer.AwsS3V4Signer
import software.amazon.awssdk.core.client.config.SdkAdvancedClientOption
import software.amazon.awssdk.http.SdkHttpConfigurationOption
import software.amazon.awssdk.http.crt.AwsCrtAsyncHttpClient
import software.amazon.awssdk.http.crt.internal.AwsCrtClientBuilderBase
import spock.lang.Specification

class S3ClientConfigurationTest extends Specification{
    def 'create S3 synchronous client configuration' (){
        given:
        def props = new Properties()
        def config = new AwsConfig([client: [connectionTimeout: 20000, maxConnections: 100, maxErrorRetry: 3, socketTimeout: 20000,
                                            proxyHost: 'host.com', proxyPort: 80, proxyScheme: 'https', proxyUsername: 'user', proxyPassword: 'pass',
                                            signerOverride: 'S3SignerType', userAgent: 'Agent1' ]])
        props.putAll(config.getS3LegacyProperties())
        when:
        def clientConfig = S3SyncClientConfiguration.create(props)
        then:
        def overrideConfig = clientConfig.getClientOverrideConfiguration()
        overrideConfig.advancedOption(SdkAdvancedClientOption.USER_AGENT_PREFIX).get() == 'Agent1'
        overrideConfig.advancedOption(SdkAdvancedClientOption.SIGNER).get() instanceof AwsS3V4Signer
        overrideConfig.retryStrategy().get().maxAttempts() == 4
        def httpClientbuilder = clientConfig.getHttpClientBuilder()
        httpClientbuilder.proxyConfiguration.host() == 'host.com'
        httpClientbuilder.proxyConfiguration.port() == 80
        httpClientbuilder.proxyConfiguration.scheme() == 'https'
        httpClientbuilder.proxyConfiguration.username() == 'user'
        httpClientbuilder.proxyConfiguration.password() == 'pass'
        httpClientbuilder.standardOptions.get(SdkHttpConfigurationOption.CONNECTION_TIMEOUT).toMillis()== 20000
        httpClientbuilder.standardOptions.get(SdkHttpConfigurationOption.READ_TIMEOUT).toMillis() == 20000 //socket timeout
        httpClientbuilder.standardOptions.get(SdkHttpConfigurationOption.MAX_CONNECTIONS) == 100
    }

    def 'create S3 asynchronous client configuration' (){
        given:
        def props = new Properties()
        def config = new AwsConfig([client: [connectionTimeout: 20000, maxConnections: 100, maxErrorRetry: 3, socketTimeout: 20000,
                                             proxyHost: 'host.com', proxyPort: 80, proxyScheme: 'https', proxyUsername: 'user', proxyPassword: 'pass',
                                             signerOverride: 'S3SignerType', userAgent: 'Agent1' ]])
        props.putAll(config.getS3LegacyProperties())
        when:
        def clientConfig = S3AsyncClientConfiguration.create(props)
        then:
        def overrideConfig = clientConfig.getClientOverrideConfiguration()
        overrideConfig.advancedOption(SdkAdvancedClientOption.USER_AGENT_PREFIX).get() == 'Agent1'
        overrideConfig.advancedOption(SdkAdvancedClientOption.SIGNER).get() instanceof AwsS3V4Signer
        overrideConfig.retryStrategy().get().maxAttempts() == 4
        def httpClientbuilder = clientConfig.getHttpClientBuilder() as AwsCrtClientBuilderBase<AwsCrtAsyncHttpClient.Builder>
        httpClientbuilder.proxyConfiguration.host() == 'host.com'
        httpClientbuilder.proxyConfiguration.port() == 80
        httpClientbuilder.proxyConfiguration.scheme() == 'https'
        httpClientbuilder.proxyConfiguration.username() == 'user'
        httpClientbuilder.proxyConfiguration.password() == 'pass'
        httpClientbuilder.getAttributeMap().get(SdkHttpConfigurationOption.CONNECTION_TIMEOUT).toMillis()== 20000
        httpClientbuilder.getAttributeMap().get(SdkHttpConfigurationOption.READ_TIMEOUT) == null //socket timeout not supported in async client
        httpClientbuilder.getAttributeMap().get(SdkHttpConfigurationOption.MAX_CONNECTIONS) == 100
    }


}
