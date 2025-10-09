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
 */

package nextflow.cloud.aws.scm

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.exception.AbortOperationException
import nextflow.scm.ProviderConfig
import software.amazon.awssdk.auth.credentials.AwsBasicCredentials
import software.amazon.awssdk.auth.credentials.AwsCredentialsProvider
import software.amazon.awssdk.auth.credentials.DefaultCredentialsProvider
import software.amazon.awssdk.auth.credentials.StaticCredentialsProvider
import software.amazon.awssdk.regions.Region

/**
 * Implements a provider config for git-remote-s3 repositories
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class S3ProviderConfig extends ProviderConfig {

    private Region region = Region.US_EAST_1

    private AwsCredentialsProvider awsCredentialsProvider = DefaultCredentialsProvider.builder().build()

    S3ProviderConfig(String name, Map values) {
        super(name, [ server: "s3://$name"] + values)
        setDefaultsFromAwsConfig()
        // Override with scm repo attributes
        setValuesFromMap(values)
    }

    S3ProviderConfig(String name){
        super(name,[ platform: 's3', server: "s3://$name"])
        setDefaultsFromAwsConfig()
    }

    private void setDefaultsFromAwsConfig() {
        final config = Global.session?.config?.aws as Map
        if( config ) {
            setValuesFromMap(config)
        }
    }
    private void setValuesFromMap(Map values){
        if( values.region ) {
            region = Region.of(values.region as String)
        }
        if( values.accessKey && values.secretKey ){
            awsCredentialsProvider = StaticCredentialsProvider.create(
                AwsBasicCredentials.builder()
                    .accessKeyId(values.accessKey as String)
                    .secretAccessKey(values.secretKey as String)
                    .build())
        }
    }

    Region getRegion(){
        this.region
    }

    AwsCredentialsProvider getAwsCredentialsProvider(){
        this.awsCredentialsProvider
    }

    @Override
    protected String resolveProjectName(String path){
        log.debug ("Resolving project name from $path. returning ")
        if (!server.startsWith('s3://'))
            new AbortOperationException("S3 project server doesn't start with s3://")
        return "${server.substring('s3://'.size())}/$path"
    }

}
