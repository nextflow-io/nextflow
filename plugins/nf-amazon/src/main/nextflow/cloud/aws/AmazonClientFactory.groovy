/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.cloud.aws

import java.nio.file.Files
import java.nio.file.Paths

import com.amazonaws.AmazonClientException
import com.amazonaws.ClientConfiguration
import com.amazonaws.auth.AWSCredentialsProvider
import com.amazonaws.auth.AWSStaticCredentialsProvider
import com.amazonaws.auth.BasicAWSCredentials
import com.amazonaws.auth.DefaultAWSCredentialsProviderChain
import com.amazonaws.auth.profile.ProfileCredentialsProvider
import com.amazonaws.client.builder.AwsClientBuilder.EndpointConfiguration
import com.amazonaws.regions.InstanceMetadataRegionProvider
import com.amazonaws.regions.Region
import com.amazonaws.regions.RegionUtils
import com.amazonaws.services.batch.AWSBatch
import com.amazonaws.services.batch.AWSBatchClient
import com.amazonaws.services.batch.AWSBatchClientBuilder
import com.amazonaws.services.ec2.AmazonEC2
import com.amazonaws.services.ec2.AmazonEC2Client
import com.amazonaws.services.ec2.AmazonEC2ClientBuilder
import com.amazonaws.services.ecs.AmazonECS
import com.amazonaws.services.ecs.AmazonECSClientBuilder
import com.amazonaws.services.logs.AWSLogs
import com.amazonaws.services.logs.AWSLogsAsyncClientBuilder
import com.amazonaws.services.s3.AmazonS3
import com.amazonaws.services.s3.AmazonS3ClientBuilder
import com.amazonaws.services.securitytoken.AWSSecurityTokenServiceClientBuilder
import com.amazonaws.services.securitytoken.model.GetCallerIdentityRequest
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.cloud.aws.config.AwsConfig
import nextflow.exception.AbortOperationException
/**
 * Implement a factory class for AWS client objects
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AmazonClientFactory {

    private AwsConfig config

    /**
     * The AWS access key credentials (optional)
     */
    private String accessKey

    /**
     * The AWS secret key credentials (optional)
     */
    private String secretKey

    /**
     * The AWS region eg. {@code eu-west-1}. If it's not specified the current region is retrieved from
     * the EC2 instance metadata
     */
    private String region

    private String profile

    /**
     * Initialise the Amazon cloud driver with default (empty) parameters
     */
    AmazonClientFactory() {
        this(new AwsConfig(Collections.emptyMap()))
    }

    AmazonClientFactory(AwsConfig config, String region=null) {
        this.config = config

        if( config.accessKey && config.secretKey ) {
            this.accessKey = config.accessKey
            this.secretKey = config.secretKey
        }

        // -- the required profile, if any
        this.profile = config.profile

        // -- check some credentials exists
        if( noCredentialsExists() )
            throw new AbortOperationException("Missing AWS security credentials -- Provide access/security keys pair or define an IAM instance profile (suggested)")

        // -- get the aws default region
        this.region = region ?: config.region ?: fetchRegion()
        if( !this.region )
            throw new AbortOperationException('Missing AWS region -- Make sure to define in your system environment the variable `AWS_DEFAULT_REGION`')
    }

    protected boolean noCredentialsExists() {
        if( this.accessKey && this.secretKey )
            return false
        if( this.profile )
            return false
        if( fetchIamRole() )
            return false
        if( SysEnv.get('AWS_ACCESS_KEY_ID') && SysEnv.get('AWS_SECRET_ACCESS_KEY') )
            return false
        final awsConfig = Paths.get(System.getProperty("user.home")).resolve(".aws/config")
        if(Files.exists(awsConfig))
            return false
        return true
    }
    /**
     * Retrieve the current IAM role eventually define for a EC2 instance.
     * See http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/iam-roles-for-amazon-ec2.html#instance-metadata-security-credentials
     *
     * @return
     *      The IAM role name associated to this instance or {@code null} if no role is defined or
     *      it's not a EC2 instance
     */
    private String fetchIamRole() {
        try {
            def stsClient = AWSSecurityTokenServiceClientBuilder.defaultClient();
            return stsClient.getCallerIdentity(new GetCallerIdentityRequest()).getArn()
        }
        catch( AmazonClientException e ) {
            log.trace "Unable to fetch IAM credentials -- Cause: ${e.message}"
            return null
        }
    }

    /**
     * Retrieve the AWS region from the EC2 instance metadata.
     * See http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-metadata.html
     *
     * @return
     *      The AWS region of the current EC2 instance eg. {@code eu-west-1} or
     *      {@code null} if it's not an EC2 instance.
     */
    private String fetchRegion() {
        try {
            return new InstanceMetadataRegionProvider().getRegion()
        }
        catch (AmazonClientException e) {
            log.debug "Cannot fetch AWS region", e
            return null
        }
    }

    /**
     * Helper method to map a region string to a {@link Region} object.
     *
     * @param region An AWS region string identifier eg. {@code eu-west-1}
     * @return A {@link Region} corresponding to the specified region string
     */
    private Region getRegionObj(String region) {
        final result = RegionUtils.getRegion(region)
        if( !result )
            throw new IllegalArgumentException("Not a valid AWS region name: $region");
        return result
    }

    /**
     * Gets or lazily creates an {@link AmazonEC2Client} instance given the current
     * configuration parameter
     *
     * @return
     *      An {@link AmazonEC2Client} instance
     */
    synchronized AmazonEC2 getEc2Client() {

        final builder = AmazonEC2ClientBuilder
                .standard()
                .withRegion(region)

        final credentials = getCredentialsProvider0()
        if( credentials )
            builder.withCredentials(credentials)

        return builder.build()
    }

    /**
     * Gets or lazily creates an {@link AWSBatchClient} instance given the current
     * configuration parameter
     *
     * @return
     *      An {@link AWSBatchClient} instance
     */
    @Memoized
    AWSBatch getBatchClient() {
        final builder = AWSBatchClientBuilder
                .standard()
                .withRegion(region)

        final credentials = getCredentialsProvider0()
        if( credentials )
            builder.withCredentials(credentials)

        return builder.build()
    }

    @Memoized
    AmazonECS getEcsClient() {

        final builder = AmazonECSClientBuilder
                .standard()
                .withRegion(region)

        final credentials = getCredentialsProvider0()
        if( credentials )
            builder.withCredentials(credentials)

        return builder.build()
    }

    @Memoized
    AWSLogs getLogsClient() {

        final builder = AWSLogsAsyncClientBuilder
                .standard()
                .withRegion(region)

        final credentials = getCredentialsProvider0()
        if( credentials )
            builder.withCredentials(credentials)

        return builder.build()
    }

    AmazonS3 getS3Client(ClientConfiguration clientConfig=null, boolean global=false) {
        final builder = AmazonS3ClientBuilder
                .standard()
                .withPathStyleAccessEnabled(config.s3Config.pathStyleAccess)
                .withForceGlobalBucketAccessEnabled(global)

        final endpoint = config.s3Config.endpoint
        if( endpoint )
            builder.withEndpointConfiguration(new EndpointConfiguration(endpoint, region))
        else
            builder.withRegion(region)

        final credentials = getCredentialsProvider0()
        if( credentials )
            builder.withCredentials(credentials)

        if( clientConfig )
            builder.withClientConfiguration(clientConfig)

        return builder.build()
    }


    protected AWSCredentialsProvider getCredentialsProvider0() {
        if( accessKey && secretKey ) {
            final creds = new BasicAWSCredentials(accessKey, secretKey)
            return new AWSStaticCredentialsProvider(creds)
        }

        if( profile ) {
            return new ProfileCredentialsProvider(profile)
        }

        return new DefaultAWSCredentialsProviderChain()
    }

}
