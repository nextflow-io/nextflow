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
import software.amazon.awssdk.services.s3.multipart.MultipartConfiguration;


import java.time.Duration;
import java.util.Properties;

/**
 * Class to convert Amazon properties in S3 asynchronous client configuration
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class S3AsyncClientConfiguration extends S3ClientConfiguration{

    private S3CrtHttpConfiguration.Builder crtHttpConfiguration;
    private MultipartConfiguration.Builder multiPartBuilder;
    private S3CrtRetryConfiguration crtRetryConfiguration;
    private Integer maxConcurrency ;
    private Double targetThroughputInGbps;
    private Long maxNativeMemoryInBytes;

    private S3CrtHttpConfiguration.Builder crtHttpConfiguration(){
        if( this.crtHttpConfiguration == null)
            this.crtHttpConfiguration = S3CrtHttpConfiguration.builder();
        return this.crtHttpConfiguration;
    }

    private MultipartConfiguration.Builder multipartBuilder(){
        if( this.multiPartBuilder == null)
            this.multiPartBuilder = MultipartConfiguration.builder();
        return this.multiPartBuilder;
    }

    public S3CrtHttpConfiguration getCrtHttpConfiguration(){
        if ( this.crtHttpConfiguration == null )
            return null;
        return this.crtHttpConfiguration.build();
    }

    public MultipartConfiguration getMultipartConfiguration(){
        if( this.multiPartBuilder == null )
            return null;
        return this.multiPartBuilder.build();
    }

    private S3AsyncClientConfiguration(){
        super();
    }

    public S3CrtRetryConfiguration getCrtRetryConfiguration(){
        return this.crtRetryConfiguration;
    }

    public Integer getMaxConcurrency(){
        return this.maxConcurrency;
    }

    public Double getTargetThroughputInGbps(){
        return this.targetThroughputInGbps;
    }

    public Long getMaxNativeMemoryInBytes(){
        return this.maxNativeMemoryInBytes;
    }

    private void setAsyncConfiguration(Properties props){

        if( props.containsKey("max_error_retry")) {
            log.trace("AWS client config - max_error_retry: {}", props.getProperty("max_error_retry"));
            this.crtRetryConfiguration = S3CrtRetryConfiguration.builder().numRetries(Integer.parseInt(props.getProperty("max_error_retry"))).build();
        }

        if( props.containsKey("max_concurrency")) {
            log.trace("AWS client config - max_concurrency: {}", props.getProperty("max_concurrency"));
            this.maxConcurrency = Integer.parseInt(props.getProperty("max_concurrency"));
        }

        if( props.containsKey("target_throughput_in_gbps")) {
            log.trace("AWS client config - target_throughput_in_gbps: {}", props.getProperty("target_throughput_in_gbps"));
            this.targetThroughputInGbps = Double.parseDouble(props.getProperty("target_throughput_in_gbps"));
        }

        if( props.containsKey("max_native_memory")) {
            log.trace("AWS client config - max_native_memory: {}", props.getProperty("max_native_memory"));
            this.maxNativeMemoryInBytes = Long.parseLong(props.getProperty("max_native_memory"));
        }

        if( props.containsKey("minimum_part_size")) {
            log.trace("AWS client config - minimum_part_size: {}", props.getProperty("minimum_part_size"));
            multipartBuilder().minimumPartSizeInBytes(Long.parseLong(props.getProperty("minimum_part_size")));
        }

        if( props.containsKey("multipart_threshold")) {
            log.trace("AWS client config - multipart_threshold: {}", props.getProperty("multipart_threshold"));
            multipartBuilder().thresholdInBytes(Long.parseLong(props.getProperty("multipart_threshold")));
        }

        if( props.containsKey("connection_timeout") ) {
            log.trace("AWS client config - connection_timeout: {}", props.getProperty("connection_timeout"));
            crtHttpConfiguration().connectionTimeout(Duration.ofMillis(Long.parseLong(props.getProperty("connection_timeout"))));
        }

        if( props.containsKey("socket_timeout")) {
            log.warn("AWS client config - 'socket_timeout' doesn't exist in AWS SDK V2 Async Client");
        }

        if( props.containsKey("proxy_host")) {
            final String host = props.getProperty("proxy_host");
            final S3CrtProxyConfiguration.Builder crtProxyConfig = S3CrtProxyConfiguration.builder();
            log.trace("AWS client config - proxy host {}", host);
            crtProxyConfig.host(host);
            if (props.containsKey("proxy_port")) {
                crtProxyConfig.port(Integer.parseInt(props.getProperty("proxy_port")));
            }
            if (props.containsKey("proxy_username")) {
                crtProxyConfig.username(props.getProperty("proxy_username"));
            }
            if (props.containsKey("proxy_password")) {
                crtProxyConfig.password(props.getProperty("proxy_password"));
            }
            if (props.containsKey("proxy_scheme")) {
                crtProxyConfig.scheme(props.getProperty("proxy_scheme"));
            }
            if (props.containsKey("proxy_domain")) {
                log.warn("AWS client config 'proxy_domain' doesn't exist in AWS SDK V2 Async Client");
            }
            if (props.containsKey("proxy_workstation")) {
                log.warn("AWS client config 'proxy_workstation' doesn't exist in AWS SDK V2 Async Client");
            }
            crtHttpConfiguration().proxyConfiguration(crtProxyConfig.build());
        }
    }

    public static S3AsyncClientConfiguration create(Properties props) {
		S3AsyncClientConfiguration config = new S3AsyncClientConfiguration();
		if( props != null ){
            config.setClientOverrideConfiguration(props);
            config.setAsyncConfiguration(props);
        }
		return config;
	}
}

