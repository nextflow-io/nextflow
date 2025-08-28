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

import software.amazon.awssdk.http.SdkHttpClient;
import software.amazon.awssdk.http.apache.ApacheHttpClient;
import software.amazon.awssdk.http.apache.ProxyConfiguration;

import java.net.URI;
import java.net.URISyntaxException;
import java.time.Duration;
import java.util.Properties;

/**
 * Class to convert Amazon properties in S3 synchronous client configuration
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class S3SyncClientConfiguration extends S3ClientConfiguration{

    private int maxConnections = 50; //Default for Sync clients
    private ApacheHttpClient.Builder httpClientBuilder;

    private ApacheHttpClient.Builder httpClientBuilder(){
        if( this.httpClientBuilder == null)
            this.httpClientBuilder = ApacheHttpClient.builder();
        return this.httpClientBuilder;
    }

    public int getMaxConnections() {
        return maxConnections;
    }

    public SdkHttpClient.Builder getHttpClientBuilder(){
        if ( this.httpClientBuilder == null )
            return null;
        return this.httpClientBuilder;
    }

    private S3SyncClientConfiguration(){
        super();
    }

    private void setClientHttpBuilder(Properties props) {
        if( props.containsKey("connection_timeout") ) {
            log.trace("AWS client config - connection_timeout: {}", props.getProperty("connection_timeout"));
            httpClientBuilder().connectionTimeout(Duration.ofMillis(Long.parseLong(props.getProperty("connection_timeout"))));
        }

        if( props.containsKey("max_connections")) {
            log.trace("AWS client config - max_connections: {}", props.getProperty("max_connections"));
            this.maxConnections = Integer.parseInt(props.getProperty("max_connections"));
            httpClientBuilder().maxConnections(this.maxConnections);
        }

        if( props.containsKey("socket_timeout")) {
            log.trace("AWS client config - socket_timeout: {}", props.getProperty("socket_timeout"));
            httpClientBuilder().socketTimeout(Duration.ofMillis(Long.parseLong(props.getProperty("socket_timeout"))));
        }

        try {
            if( props.containsKey("proxy_host")) {
                final String host = props.getProperty("proxy_host");
                final int port = Integer.parseInt(props.getProperty("proxy_port", "-1"));
                final String scheme = props.getProperty("proxy_scheme", "http");
                final ProxyConfiguration.Builder proxyConfig = ProxyConfiguration.builder();
                log.trace("AWS client config - proxy {}://{}:{}", scheme, host, port);
                proxyConfig.endpoint(new URI(scheme, null, host, port, null, null, null));

                if (props.containsKey("proxy_username")) {
                    proxyConfig.username(props.getProperty("proxy_username"));
                }
                if (props.containsKey("proxy_password")) {
                    proxyConfig.password(props.getProperty("proxy_password"));
                }

                if (props.containsKey("proxy_domain")) {
                    proxyConfig.ntlmDomain(props.getProperty("proxy_domain"));
                }
                if (props.containsKey("proxy_workstation")) {
                    proxyConfig.ntlmWorkstation(props.getProperty("proxy_workstation"));
                }

                httpClientBuilder().proxyConfiguration(proxyConfig.build());
            }
        } catch (URISyntaxException e){
            log.warn("Exception creating AWS client config - proxy URI", e);
        }
    }

    public static S3SyncClientConfiguration create(Properties props) {
		S3SyncClientConfiguration config = new S3SyncClientConfiguration();

		if( props != null ) {
            config.setClientOverrideConfiguration(props);
            config.setClientHttpBuilder(props);
        }

		return config;
	}


}

