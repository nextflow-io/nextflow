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

import nextflow.cloud.aws.config.AwsBucketConfig
import nextflow.cloud.aws.nio.util.S3AsyncClientConfiguration
import nextflow.cloud.aws.nio.util.S3SyncClientConfiguration
import software.amazon.awssdk.auth.credentials.AnonymousCredentialsProvider
import software.amazon.awssdk.auth.credentials.AwsCredentialsProvider
import software.amazon.awssdk.auth.credentials.DefaultCredentialsProvider
import software.amazon.awssdk.auth.credentials.StaticCredentialsProvider
import software.amazon.awssdk.auth.credentials.AwsBasicCredentials
import software.amazon.awssdk.auth.credentials.ProfileCredentialsProvider
import software.amazon.awssdk.core.client.config.ClientOverrideConfiguration
import software.amazon.awssdk.core.exception.SdkClientException
import software.amazon.awssdk.http.SdkHttpClient
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.regions.providers.InstanceProfileRegionProvider
import software.amazon.awssdk.services.batch.BatchClient
import software.amazon.awssdk.services.cloudwatchlogs.CloudWatchLogsClient
import software.amazon.awssdk.services.ec2.Ec2Client
import software.amazon.awssdk.services.ecs.EcsClient
import software.amazon.awssdk.services.s3.S3AsyncClient
import software.amazon.awssdk.services.s3.S3Client
import software.amazon.awssdk.services.s3.S3Configuration
import software.amazon.awssdk.services.s3.S3CrtAsyncClientBuilder
import software.amazon.awssdk.services.s3.multipart.MultipartConfiguration
import software.amazon.awssdk.services.sts.StsClient
import software.amazon.awssdk.services.sts.model.GetCallerIdentityRequest
import software.amazon.awssdk.services.sts.model.StsException
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
        try {
            final stsClient = StsClient.create()
            return stsClient.getCallerIdentity(GetCallerIdentityRequest.builder().build() as GetCallerIdentityRequest).arn();
        }
        catch (StsException e) {
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
        }
        catch (SdkClientException e) {
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
        return Ec2Client.builder()
            .region(getRegionObj(region))
            .credentialsProvider(getCredentialsProvider0())
            .build()
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
        return BatchClient.builder()
            .region(getRegionObj(region))
            .credentialsProvider(getCredentialsProvider0())
            .build()
    }

    @Memoized
    EcsClient getEcsClient() {
        return EcsClient.builder()
            .region(getRegionObj(region))
            .credentialsProvider(getCredentialsProvider0())
            .build()
    }

    @Memoized
    CloudWatchLogsClient getLogsClient() {
        return CloudWatchLogsClient.builder().region(getRegionObj(region)).credentialsProvider(getCredentialsProvider0()).build()
    }

    S3Client getS3Client(S3SyncClientConfiguration s3ClientConfig, String bucketName, boolean global = false) {
        final SdkHttpClient.Builder httpClientBuilder = s3ClientConfig.getHttpClientBuilder()
        final ClientOverrideConfiguration overrideConfiguration = s3ClientConfig.getClientOverrideConfiguration()
        final bucketConfig = config.getBucketConfig(bucketName)
        final builder = S3Client.builder()
            .crossRegionAccessEnabled(global)
            .credentialsProvider(getS3CredentialsProvider(bucketConfig))
            .serviceConfiguration(S3Configuration.builder()
                .pathStyleAccessEnabled(bucketConfig.s3PathStyleAccess)
                .multiRegionEnabled(global)
                .build())

        if( bucketConfig.endpoint )
            builder.endpointOverride(URI.create(bucketConfig.endpoint))

        // AWS SDK v2 region must be always set, even when endpoint is overridden
        builder.region(getRegionObj(region))

        if( httpClientBuilder != null )
            builder.httpClientBuilder(httpClientBuilder)

        if( overrideConfiguration != null )
            builder.overrideConfiguration(overrideConfiguration)

        return builder.build()
    }

    S3AsyncClient getS3AsyncClient(S3AsyncClientConfiguration s3ClientConfig, String bucketName, boolean global = false) {
        final bucketConfig = config.getBucketConfig(bucketName)
        def builder = S3AsyncClient.crtBuilder()
            .crossRegionAccessEnabled(global)
            .credentialsProvider(getS3CredentialsProvider(bucketConfig))
            .forcePathStyle(bucketConfig.pathStyleAccess)
            .region(getRegionObj(region))
        if( bucketConfig.endpoint )
            builder.endpointOverride(URI.create(bucketConfig.endpoint))

        final retryConfiguration = s3ClientConfig.getCrtRetryConfiguration()
        if( retryConfiguration != null )
            builder.retryConfiguration(retryConfiguration)

        final httpConfiguration = s3ClientConfig.getCrtHttpConfiguration()
        if( httpConfiguration != null )
            builder.httpConfiguration(httpConfiguration)

        final multipartConfig = s3ClientConfig.getMultipartConfiguration()
        if( multipartConfig != null )
            setMultipartConfiguration(multipartConfig, builder)

        final throughput = s3ClientConfig.getTargetThroughputInGbps()
        if( throughput != null )
            builder.targetThroughputInGbps(throughput)

        final nativeMemory = s3ClientConfig.getMaxNativeMemoryInBytes()
        if (nativeMemory != null )
            builder.maxNativeMemoryLimitInBytes(nativeMemory)

        final maxConcurrency = s3ClientConfig.getMaxConcurrency()
        if( maxConcurrency != null )
            builder.maxConcurrency(maxConcurrency)

        return builder.build()
    }
    /**
     * Returns an AwsCredentialsProvider for S3 clients.
     *
     * This method wraps the same AWS credentials used for other clients, but ensures proper handling of anonymous S3 access.
     * If the 'anonymous' flag is set in Nextflow's AWS S3 configuration, or if no credentials are resolved by other providers,
     * an AnonymousCredentialsProvider instance is returned.
     *
     * Prior to AWS SDK v2, the S3CredentialsProvider automatically managed fallback to anonymous access when no credentials were found.
     * However, due to a limitation in the AWS SDK v2 CRT Async S3 client (see https://github.com/aws/aws-sdk-java-v2/issues/5810),
     * anonymous credentials only work when explicitly configured via AnonymousCredentialsProvider.
     * Custom credential providers or provider chains that resolve to anonymous credentials are not handled correctly by the CRT client.
     *
     * To work around this, this method explicitly checks whether credentials can be resolved.
     * If no credentials are found, it returns an AnonymousCredentialsProvider; otherwise, it returns the resolved provider.
     *
     * @return an AwsCredentialsProvider instance, falling back to anonymous if needed.
     */
    private AwsCredentialsProvider getS3CredentialsProvider(AwsBucketConfig configBucket) {
        if ( configBucket?.anonymous || config.s3Config.anonymous)
            return AnonymousCredentialsProvider.create()
        def provider = getCredentialsProvider0()
        try {
            provider.resolveCredentials()
        } catch (Exception e) {
            log.debug("No AWS credentials available - falling back to anonymous access")
            return AnonymousCredentialsProvider.create()
        }
        return provider
    }

    private void setMultipartConfiguration(MultipartConfiguration multipartConfig, S3CrtAsyncClientBuilder builder) {
        if( multipartConfig.minimumPartSizeInBytes() != null )
            builder.minimumPartSizeInBytes(multipartConfig.minimumPartSizeInBytes())
        if( multipartConfig.thresholdInBytes() != null )
            builder.thresholdInBytes(multipartConfig.thresholdInBytes())
    }

    protected AwsCredentialsProvider getCredentialsProvider0() {
        if( accessKey && secretKey ) {
            return StaticCredentialsProvider.create(AwsBasicCredentials.create(accessKey, secretKey))
        }

        if( profile ) {
            return ProfileCredentialsProvider.builder()
                    .profileName(profile)
                    .build()
        }

        return DefaultCredentialsProvider.create()
    }

}
