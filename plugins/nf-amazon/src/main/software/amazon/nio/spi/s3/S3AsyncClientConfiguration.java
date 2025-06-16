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
package software.amazon.nio.spi.s3;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import software.amazon.awssdk.auth.signer.Aws4Signer;
import software.amazon.awssdk.auth.signer.AwsS3V4Signer;
import software.amazon.awssdk.core.client.config.ClientOverrideConfiguration;
import software.amazon.awssdk.core.client.config.SdkAdvancedClientOption;
import software.amazon.awssdk.core.signer.Signer;
import software.amazon.awssdk.services.s3.crt.S3CrtHttpConfiguration;
import software.amazon.awssdk.services.s3.crt.S3CrtProxyConfiguration;
import software.amazon.awssdk.services.s3.crt.S3CrtRetryConfiguration;

import java.time.Duration;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

/**
 * Class to convert Amazon properties in S3 asynchronous client configuration
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class S3AsyncClientConfiguration {

    private static final Logger log = LoggerFactory.getLogger(S3AsyncClientConfiguration.class);

    private S3CrtRetryConfiguration.Builder retryConfiguration;
    private S3CrtHttpConfiguration.Builder httpConfiguration;
    private ClientOverrideConfiguration.Builder cocBuilder;
    private NextflowS3ClientOpenOptions openOptions;
    private Integer maxConcurrency;

    private S3CrtRetryConfiguration.Builder retryConfiguration(){
        if( this.retryConfiguration == null )
            this.retryConfiguration = S3CrtRetryConfiguration.builder();
        return this.retryConfiguration;
    }
    
    public boolean hasRetryConfiguration() { return this.retryConfiguration != null; }
    
    private S3CrtHttpConfiguration.Builder httpConfiguration(){
        if( this.httpConfiguration == null)
            this.httpConfiguration = S3CrtHttpConfiguration.builder();
        return this.httpConfiguration;
    }
    
    public boolean hasHttpConfiguration(){ return this.httpConfiguration != null; }

    private ClientOverrideConfiguration.Builder cocBuilder(){
        if( this.cocBuilder == null )
            this.cocBuilder = ClientOverrideConfiguration.builder();
        return this.cocBuilder;
    }
    
    private NextflowS3ClientOpenOptions openOptions(){
        if( this.openOptions == null )
            this.openOptions = new NextflowS3ClientOpenOptions();
        return this.openOptions;
    }
    
    public boolean hasClientOverrideConfiguration() { return this.cocBuilder != null; }
    
    
    public S3CrtRetryConfiguration getRetryConfiguration(){
        if( retryConfiguration == null )
            return null;
        return retryConfiguration.build();
    }
    
    public S3CrtHttpConfiguration getHttpConfiguration(){
        if( httpConfiguration == null )
            return null;
        return httpConfiguration.build();
    }
    
    public ClientOverrideConfiguration getClientOverrideConfiguration(){
        if ( this.cocBuilder == null )
            return null;
        return this.cocBuilder.build();
    }
    
    public Integer getMaxConcurrency() { return this.maxConcurrency; }

    public void setMaxConcurrency(Integer maxConcurrency) { 
        if( maxConcurrency == null )
            return;
        this.maxConcurrency = maxConcurrency;
    }
    
    public boolean hasMaxConcurrency() { return this.maxConcurrency != null; }
    
    private S3AsyncClientConfiguration(){}
    
    public static S3AsyncClientConfiguration create(Properties props) {
		S3AsyncClientConfiguration config = new S3AsyncClientConfiguration();

		if( props == null )
			return config;

		if( props.containsKey("connection_timeout") ) {
			log.trace("AWS client config - connection_timeout: {}", props.getProperty("connection_timeout"));
			config.httpConfiguration().connectionTimeout(Duration.ofMillis(Long.parseLong(props.getProperty("connection_timeout"))));
		}
        
        if( props.containsKey("max_connections")) {
			log.trace("AWS client config - max_connections: {}", props.getProperty("max_connections"));
			config.setMaxConcurrency(Integer.parseInt(props.getProperty("max_connections")));
		}

		if( props.containsKey("max_error_retry")) {
			log.trace("AWS client config - max_error_retry: {}", props.getProperty("max_error_retry"));
			config.retryConfiguration().numRetries(Integer.parseInt(props.getProperty("max_error_retry")));
		}

		if( props.containsKey("protocol")) {
			log.warn("AWS client config 'protocol' doesn't exist in AWS SDK V2");
		}
        if( props.containsKey("proxy_host")) {
            setProxyProperties(props, config);
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

    public List<S3OpenOption> getOpenOptions() {
        final List<S3OpenOption> list = new ArrayList<S3OpenOption>();
        if( this.cocBuilder != null )
            openOptions().setClientOverride(this.getClientOverrideConfiguration());
        if( this.openOptions != null ) {
            list.add(this.openOptions.copy());
        }
        return list;
    }

    private static void setOpenOptions(Properties props, S3AsyncClientConfiguration config) {
        if( props.containsKey("requester_pays_enabled") ) {
            log.trace("AWS client config - requester_pays_enabled : {}", props.getProperty("requester_pays_enabled"));
            config.openOptions().setRequesterPays(props.getProperty("requester_pays_enabled"));
        }
        String aclProp = getProp(props, "s_3_acl", "s3_acl", "s3Acl");
        if( aclProp != null && !aclProp.isBlank()) {
            log.trace("AWS client config - acl : {}", aclProp);
            config.openOptions().setCannedAcl(aclProp);
        }
        if( props.containsKey("storage_encryption") ){
            log.trace("AWS client config - storage_encryption : {}", props.getProperty("storage_encryption"));
            config.openOptions().setStorageEncryption(props.getProperty("storage_encryption"));
        }
        if( props.containsKey("storage_kms_key_id") ) {
            log.trace("AWS client config - storage_kms_key_id : {}", props.getProperty("storage_kms_key_id"));
            config.openOptions().setKmsKeyId(props.getProperty("storage_kms_key_id"));
        }

    }

    private static String getProp(Properties props, String... keys) {
        for( String k : keys ) {
            if( props.containsKey(k) ) {
                return props.getProperty(k);
            }
        }
        return null;
    }

    private static void setProxyProperties(Properties props, S3AsyncClientConfiguration config) {
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

        if (props.containsKey("proxy_domain")) {
            log.warn("AWS client config 'proxy_domain' doesn't exist in AWS SDK V2 Async Client");
        }
        if (props.containsKey("proxy_workstation")) {
            log.warn("AWS client config 'proxy_workstation' doesn't exist in AWS SDK V2 Async Client");
        }

        config.httpConfiguration().proxyConfiguration(proxyConfig.build());
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

