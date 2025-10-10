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
import nextflow.cloud.aws.AwsClientFactory
import nextflow.cloud.aws.config.AwsConfig
import nextflow.exception.AbortOperationException
import nextflow.scm.ProviderConfig
import software.amazon.awssdk.auth.credentials.AwsBasicCredentials
import software.amazon.awssdk.auth.credentials.AwsCredentialsProvider
import software.amazon.awssdk.auth.credentials.DefaultCredentialsProvider
import software.amazon.awssdk.auth.credentials.ProfileCredentialsProvider
import software.amazon.awssdk.auth.credentials.StaticCredentialsProvider
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.regions.providers.DefaultAwsRegionProviderChain

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
        setValues(values)
    }

    S3ProviderConfig(String name){
        super(name,[ platform: 's3', server: "s3://$name"])
        setValues()
    }

    private void setValues(Map values = Map.of()) {
        //Get sessions config if exists
        def session = Global.session?.config?.aws as Map ?: Map.of()

        //Merge with scm values and convert to AwsConfg to unify SysEnv fallback and profile management
        final config = new AwsConfig(session + values)
        if( config.region ) {
            region = Region.of(config.region)
        }
        if( config.accessKey && config.secretKey ){
            awsCredentialsProvider = StaticCredentialsProvider.create(
                AwsBasicCredentials.builder()
                    .accessKeyId(config.accessKey as String)
                    .secretAccessKey(config.secretKey as String)
                    .build())
        } else if( config.profile ){
            // Get credentials from profile
            awsCredentialsProvider = ProfileCredentialsProvider.builder().profileName(config.profile).build()
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
