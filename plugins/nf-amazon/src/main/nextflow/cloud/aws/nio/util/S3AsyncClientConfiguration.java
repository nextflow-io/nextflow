/*
 * Copyright 2020-2025, Seqera Labs
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
package nextflow.cloud.aws.nio.util;

import nextflow.cloud.aws.nio.S3Client;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import software.amazon.awssdk.auth.signer.Aws4Signer;
import software.amazon.awssdk.auth.signer.AwsS3V4Signer;
import software.amazon.awssdk.core.client.config.ClientOverrideConfiguration;
import software.amazon.awssdk.core.client.config.SdkAdvancedClientOption;
import software.amazon.awssdk.core.signer.Signer;
import software.amazon.awssdk.http.async.SdkAsyncHttpClient;
import software.amazon.awssdk.http.crt.AwsCrtAsyncHttpClient;
import software.amazon.awssdk.http.crt.ProxyConfiguration;

import java.time.Duration;
import java.util.Properties;

/**
 * Class to convert Amazon properties in S3 asynchronous client configuration
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class S3AsyncClientConfiguration {

    private static final Logger log = LoggerFactory.getLogger(S3Client.class);

    private ClientOverrideConfiguration.Builder cocBuilder;
    private AwsCrtAsyncHttpClient.Builder httpClientBuilder;

    private ClientOverrideConfiguration.Builder cocBuilder(){
        if( this.cocBuilder == null )
            this.cocBuilder = ClientOverrideConfiguration.builder();
        return this.cocBuilder;
    }

    private AwsCrtAsyncHttpClient.Builder httpClientBuilder(){
        if( this.httpClientBuilder == null)
            this.httpClientBuilder = AwsCrtAsyncHttpClient.builder();
        return this.httpClientBuilder;
    }

    public ClientOverrideConfiguration getClientOverrideConfiguration(){
        if( cocBuilder == null )
            return null;
        return cocBuilder.build();
    }
     public SdkAsyncHttpClient getHttpClient(){
        if ( this.httpClientBuilder == null )
            return null;
        return this.httpClientBuilder.build();
     }

    private S3AsyncClientConfiguration(){}


    public static S3AsyncClientConfiguration create(Properties props) {
		S3AsyncClientConfiguration config = new S3AsyncClientConfiguration();

		if( props == null )
			return config;

		if( props.containsKey("connection_timeout") ) {
			log.trace("AWS client config - connection_timeout: {}", props.getProperty("connection_timeout"));
			config.httpClientBuilder().connectionTimeout(Duration.ofMillis(Long.parseLong(props.getProperty("connection_timeout"))));
		}

		if( props.containsKey("max_connections")) {
			log.trace("AWS client config - max_connections: {}", props.getProperty("max_connections"));
			config.httpClientBuilder().maxConcurrency(Integer.parseInt(props.getProperty("max_connections")));
		}

		if( props.containsKey("max_error_retry")) {
			log.trace("AWS client config - max_error_retry: {}", props.getProperty("max_error_retry"));
			config.cocBuilder().retryStrategy(b -> b.maxAttempts(Integer.parseInt(props.getProperty("max_error_retry"))));
		}

		if( props.containsKey("protocol")) {
			log.warn("AWS client config 'protocol' doesn't exist in AWS SDK V2");
		}
        if( props.containsKey("proxy_host")) {
            final String host = props.getProperty("proxy_host");
            final ProxyConfiguration.Builder proxyConfig = ProxyConfiguration.builder();
            log.trace("AWS client config - proxy host {}", host);
            proxyConfig.host(host);
            if (props.containsKey("proxy_port")) {
                proxyConfig.port(Integer.parseInt(props.getProperty("proxy_port")));
            }
            if (props.containsKey("proxy_username")) {
                proxyConfig.username(props.getProperty("proxy_username"));
            }
            if (props.containsKey("proxy_password")) {
                proxyConfig.password(props.getProperty("proxy_password"));
            }

            if (props.containsKey("proxy_domain")) {
                log.warn("AWS client config 'proxy_domain' doesn't exist in AWS SDK V2 Async Client");
            }
            if (props.containsKey("proxy_workstation")) {
                log.warn("AWS client config 'proxy_workstation' doesn't exist in AWS SDK V2 Async Client");
            }

            config.httpClientBuilder().proxyConfiguration(proxyConfig.build());
        }

		if ( props.containsKey("signer_override")) {
			log.warn("AWS client config - 'signerOverride' is deprecated");
			config.cocBuilder().putAdvancedOption(SdkAdvancedClientOption.SIGNER, resolveSigner(props.getProperty("signer_override")));
		}

		if( props.containsKey("socket_send_buffer_size_hints") || props.containsKey("socket_recv_buffer_size_hints") ) {
			log.warn("AWS client config - 'socket_send_buffer_size_hints' and 'socket_recv_buffer_size_hints' do not exist in AWS SDK V2" );
		}

		if( props.containsKey("socket_timeout")) {
			log.warn("AWS client config - 'socket_timeout' doesn't exist in AWS SDK V2 Async Client");

		}

		if( props.containsKey("user_agent")) {
			log.trace("AWS client config - user_agent: {}", props.getProperty("user_agent"));
			config.cocBuilder().putAdvancedOption(SdkAdvancedClientOption.USER_AGENT_PREFIX, props.getProperty("user_agent"));
		}
		return config;
	}

    private static Signer resolveSigner(String signerOverride) {
        switch (signerOverride) {
            case "AWSS3V4SignerType":
            case "S3SignerType":
                return AwsS3V4Signer.create();
        case "AWS4SignerType":
            return Aws4Signer.create();
        default:
            throw new IllegalArgumentException("Unsupported signer: " + signerOverride);
    }
}
}

