/*
 * Copyright 2013-2024, Seqera Labs
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

import software.amazon.awssdk.auth.credentials.AnonymousCredentialsProvider
import software.amazon.awssdk.auth.credentials.AwsCredentialsProvider
import software.amazon.awssdk.auth.credentials.DefaultCredentialsProvider
import software.amazon.awssdk.auth.credentials.StaticCredentialsProvider
import software.amazon.awssdk.auth.credentials.AwsBasicCredentials
import software.amazon.awssdk.auth.credentials.ProfileCredentialsProvider
import software.amazon.awssdk.core.exception.SdkClientException
import software.amazon.awssdk.http.SdkHttpClient
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.regions.providers.InstanceProfileRegionProvider
import software.amazon.awssdk.services.batch.BatchClient
import software.amazon.awssdk.services.cloudwatchlogs.CloudWatchLogsClient
import software.amazon.awssdk.services.ec2.Ec2Client
import software.amazon.awssdk.services.ecs.EcsClient
import software.amazon.awssdk.services.s3.S3Client
import software.amazon.awssdk.services.s3.S3Configuration
import software.amazon.awssdk.services.sts.StsClient
import software.amazon.awssdk.services.sts.model.GetCallerIdentityRequest
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.cloud.aws.config.AwsConfig
import nextflow.cloud.aws.util.S3CredentialsProvider

import nextflow.exception.AbortOperationException
import software.amazon.awssdk.services.sts.model.StsException

/**
 * Implement a factory class for AWS client objects
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsClientFactory {

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
    AwsClientFactory() {
        this(new AwsConfig(Collections.emptyMap()))
    }

    AwsClientFactory(AwsConfig config, String region=null) {
        this.config = config

        if( config.accessKey && config.secretKey ) {
            this.accessKey = config.accessKey
            this.secretKey = config.secretKey
        }

        // -- the required profile, if any
        this.profile
                = config.profile
                ?: SysEnv.get('AWS_PROFILE')
                ?: SysEnv.get('AWS_DEFAULT_PROFILE')

        // -- get the aws default region
        this.region
                = region
                ?: config.region
                ?: SysEnv.get('AWS_REGION')
                ?: SysEnv.get('AWS_DEFAULT_REGION')
                ?: fetchRegion()

        if( !this.region )
            throw new AbortOperationException('Missing AWS region -- Make sure to define in your system environment the variable `AWS_DEFAULT_REGION`')
    }

    String accessKey() { accessKey }

    String secretKey() { secretKey }

    String region() { region }

    String profile() { profile }

    /**
     * Retrieve the current IAM role eventually define for a EC2 instance.
     * See http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/iam-roles-for-amazon-ec2.html#instance-metadata-security-credentials
     *
     * @return
     *      The IAM role name associated to this instance or {@code null} if no role is defined or
     *      it's not a EC2 instance
     */
    protected String fetchIamRole() {
        try{
            StsClient stsClient = StsClient.create()
            return stsClient.getCallerIdentity(GetCallerIdentityRequest.builder().build() as GetCallerIdentityRequest).arn();
        } catch ( StsException e) {
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
            return new InstanceProfileRegionProvider().getRegion().id();
        } catch ( SdkClientException e) {
            log.debug("Cannot fetch AWS region", e);
            return null;
        }
    }

    /**
     * Helper method to map a region string to a {@link Region} object.
     *
     * @param region An AWS region string identifier eg. {@code eu-west-1}
     * @return A {@link Region} corresponding to the specified region string
     */
    private Region getRegionObj(String region) {
        final result = Region.of(region)
        if( !result )
            throw new IllegalArgumentException("Not a valid AWS region name: $region");
        return result
    }

    /**
     * Gets or lazily creates an {@link Ec2Client} instance given the current
     * configuration parameter
     *
     * @return
     *      An {@link Ec2Client} instance
     */
    synchronized Ec2Client getEc2Client() {
        return Ec2Client.builder().region(getRegionObj(region)).credentialsProvider(getCredentialsProvider0()).build()
    }

    /**
     * Gets or lazily creates an {@link BatchClient} instance given the current
     * configuration parameter
     *
     * @return
     *      An {@link BatchClient} instance
     */
    @Memoized
    BatchClient getBatchClient() {
        return BatchClient.builder().region(getRegionObj(region)).credentialsProvider(getCredentialsProvider0()).build()
    }

    @Memoized
    EcsClient getEcsClient() {
        return EcsClient.builder().region(getRegionObj(region)).credentialsProvider(getCredentialsProvider0()).build()
    }

    @Memoized
    CloudWatchLogsClient getLogsClient() {
        return CloudWatchLogsClient.builder().region(getRegionObj(region)).credentialsProvider(getCredentialsProvider0()).build()
    }

    S3Client getS3Client(SdkHttpClient httpClient = null, boolean global = false) {
        def builder = S3Client.builder()
            .credentialsProvider(config.s3Config.anonymous ? AnonymousCredentialsProvider.create() : new S3CredentialsProvider(getCredentialsProvider0()))
            .serviceConfiguration(S3Configuration.builder()
                .pathStyleAccessEnabled(config.s3Config.pathStyleAccess)
                .build())

        if (config.s3Config.endpoint) {
            builder.endpointOverride(URI.create(config.s3Config.endpoint))
        } else {
            builder.region(getRegionObj(region))
        }

        if (httpClient != null) {
            builder.httpClient(httpClient)
        }

        return builder.build()
    }

    protected AwsCredentialsProvider getCredentialsProvider0() {
        if (accessKey && secretKey) {
            return StaticCredentialsProvider.create(AwsBasicCredentials.create(accessKey, secretKey))
        }

        if (profile) {
            return ProfileCredentialsProvider.builder()
                    .profileName(profile)
                    .build()
        }

        return DefaultCredentialsProvider.create()
    }

}
