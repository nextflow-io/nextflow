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

import software.amazon.awssdk.services.s3.crt.S3CrtProxyConfiguration;
import software.amazon.awssdk.services.s3.crt.S3CrtHttpConfiguration;
import software.amazon.awssdk.services.s3.crt.S3CrtRetryConfiguration;


import java.time.Duration;
import java.util.Properties;

/**
 * Class to convert Amazon properties in S3 asynchronous client configuration
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class S3AsyncClientConfiguration extends S3ClientConfiguration{

    private S3CrtHttpConfiguration.Builder httpConfiguration;
    private S3CrtRetryConfiguration retryConfiguration;
    private Integer maxConcurrency;

    private S3CrtHttpConfiguration.Builder httpClientBuilder(){
        if( this.httpConfiguration == null)
            this.httpConfiguration = S3CrtHttpConfiguration.builder();
        return this.httpConfiguration;
    }

     public S3CrtHttpConfiguration getHttpConfiguration(){
        if ( this.httpConfiguration == null )
            return null;
        return this.httpConfiguration.build();
     }

    private S3AsyncClientConfiguration(){
        super();
    }

    private void setMaxConcurrency(Properties props){
        if( props.containsKey("max_connections")) {
            log.trace("AWS client config - max_connections: {}", props.getProperty("max_connections"));
            this.maxConcurrency = Integer.parseInt(props.getProperty("max_connections"));
        }
    }

    private void setRetryConfiguration(Properties props) {

		if( props.containsKey("max_error_retry")) {
			log.trace("AWS client config - max_error_retry: {}", props.getProperty("max_error_retry"));
			this.retryConfiguration = S3CrtRetryConfiguration.builder().numRetries(Integer.parseInt(props.getProperty("max_error_retry"))).build();
		}
    }

    public S3CrtRetryConfiguration getRetryConfiguration(){
        return this.retryConfiguration;
    }

    public Integer getMaxConcurrency(){
        return this.maxConcurrency;
    }

    private void setHttpConfiguration(Properties props){
        if( props.containsKey("connection_timeout") ) {
            log.trace("AWS client config - connection_timeout: {}", props.getProperty("connection_timeout"));
            httpClientBuilder().connectionTimeout(Duration.ofMillis(Long.parseLong(props.getProperty("connection_timeout"))));
        }

        if( props.containsKey("socket_timeout")) {
            log.warn("AWS client config - 'socket_timeout' doesn't exist in AWS SDK V2 Async Client");
        }

        if( props.containsKey("proxy_host")) {
            final String host = props.getProperty("proxy_host");
            final S3CrtProxyConfiguration.Builder proxyConfig = S3CrtProxyConfiguration.builder();
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
            if (props.containsKey("proxy_scheme")) {
                proxyConfig.scheme(props.getProperty("proxy_scheme"));
            }
            if (props.containsKey("proxy_domain")) {
                log.warn("AWS client config 'proxy_domain' doesn't exist in AWS SDK V2 Async Client");
            }
            if (props.containsKey("proxy_workstation")) {
                log.warn("AWS client config 'proxy_workstation' doesn't exist in AWS SDK V2 Async Client");
            }
            httpClientBuilder().proxyConfiguration(proxyConfig.build());
        }

    }

    public static S3AsyncClientConfiguration create(Properties props) {
		S3AsyncClientConfiguration config = new S3AsyncClientConfiguration();
		if( props != null ){
            config.setClientOverrideConfiguration(props);
            config.setHttpConfiguration(props);
            config.setRetryConfiguration(props);
        }
		return config;
	}
}

