/*
 * Copyright 2013-2026, Seqera Labs
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
 */
package nextflow.cloud.aws.nio.util;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import software.amazon.awssdk.auth.signer.Aws4Signer;
import software.amazon.awssdk.auth.signer.AwsS3V4Signer;
import software.amazon.awssdk.core.client.config.ClientOverrideConfiguration;
import software.amazon.awssdk.core.client.config.SdkAdvancedClientOption;
import software.amazon.awssdk.core.signer.Signer;
import software.amazon.awssdk.retries.StandardRetryStrategy;
import software.amazon.awssdk.utils.ClassLoaderHelper;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Properties;

/**
 * Class to convert Amazon properties in S3 client override configuration
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class S3ClientConfiguration {

    protected static final Logger log = LoggerFactory.getLogger(S3ClientConfiguration.class);

    private ClientOverrideConfiguration.Builder cocBuilder;

    private ClientOverrideConfiguration.Builder cocBuilder(){
        if( this.cocBuilder == null )
            this.cocBuilder = ClientOverrideConfiguration.builder();
        return this.cocBuilder;
    }

    public ClientOverrideConfiguration getClientOverrideConfiguration(){
        if( cocBuilder == null )
            return null;
        return cocBuilder.build();
    }

    protected S3ClientConfiguration(){}


    protected final void setClientOverrideConfiguration(Properties props) {
		if( props == null )
			return;

		if( props.containsKey("max_error_retry")) {
			log.trace("AWS client config - max_error_retry: {}", props.getProperty("max_error_retry"));
			cocBuilder().retryStrategy(StandardRetryStrategy.builder().maxAttempts((Integer.parseInt(props.getProperty("max_error_retry")) + 1 )).build());
		}

		if( props.containsKey("protocol")) {
			log.warn("AWS client config 'protocol' doesn't exist in AWS SDK V2");
		}

		if ( props.containsKey("signer_override")) {
			log.warn("AWS client config 'signerOverride' is not supported in AWS SDK V2. This option will be ignored.");

		}

		if( props.containsKey("socket_send_buffer_size_hints") || props.containsKey("socket_recv_buffer_size_hints") ) {
			log.warn("AWS client config - 'socket_send_buffer_size_hints' and 'socket_recv_buffer_size_hints' do not exist in AWS SDK V2" );
		}

		if( props.containsKey("user_agent")) {
			log.warn("AWS client config 'user_agent' is not supported in AWS SDK V2. This option will be ignored.");
		}
	}
}

