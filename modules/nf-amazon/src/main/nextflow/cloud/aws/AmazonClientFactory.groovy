/*
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

import com.amazonaws.auth.AWSCredentials
import com.amazonaws.auth.AWSStaticCredentialsProvider
import com.amazonaws.auth.BasicAWSCredentials
import com.amazonaws.auth.BasicSessionCredentials
import com.amazonaws.regions.Region
import com.amazonaws.regions.RegionUtils
import com.amazonaws.services.batch.AWSBatchClient
import com.amazonaws.services.ec2.AmazonEC2Client
import com.amazonaws.services.ecs.AmazonECS
import com.amazonaws.services.ecs.AmazonECSClientBuilder
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.exception.AbortOperationException
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AmazonClientFactory {

    /**
     * Reference to {@link AmazonEC2Client} object
     */
    private AmazonEC2Client ec2Client

    /**
     * The AWS access key credentials (optional)
     */
    private String accessKey

    /**
     * The AWS secret key credentials (optional)
     */
    private String secretKey

    /**
     * The AWS session key credentials (optional)
     */
    private String sessionToken

    /**
     * The AWS region eg. {@code eu-west-1}. If it's not specified the current region is retrieved from
     * the EC2 instance metadata
     */
    private String region

    /**
     * @return The current set AWS access key
     */
    String getAccessKey() { accessKey }

    /**
     * @return The current set AWS secret key
     */
    String getSecretKey() { secretKey }


    String getSessionToken() { sessionToken }


    /**
     * Initialise the Amazon cloud driver with default (empty) parameters
     */
    AmazonClientFactory() {
        this(Collections.emptyMap())
    }

    /**
     * Initialise the Amazon cloud driver with the specified parameters
     *
     * @param config
     *      A map holding the driver parameters:
     *      - accessKey: the access key credentials
     *      - secretKey: the secret key credentials
     *      - region: the AWS region
     */
    AmazonClientFactory(Map config) {
        // -- get the aws credentials
        List credentials
        if( config.accessKey && config.secretKey ) {
            this.accessKey = config.accessKey
            this.secretKey = config.secretKey
            if (config.sessionToken){
                this.sessionToken = config.sessionToken
            }
        }
        else if( (credentials= Global.getAwsCredentials()) ) {
            this.accessKey = credentials[0]
            this.secretKey = credentials[1]
            if (credentials.size() == 3){
                this.sessionToken = credentials[2]
            }
        }

        if( !accessKey && !fetchIamRole() )
            throw new AbortOperationException("Missing AWS security credentials -- Provide access/security keys pair or define a IAM instance profile (suggested)")

        // -- get the aws default region
        region = config.region ?: Global.getAwsRegion() ?: fetchRegion()
        if( !region )
            throw new AbortOperationException('Missing AWS region -- Make sure to define in your system environment the variable `AWS_DEFAULT_REGION`')

    }

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
            def role = getUrl('http://169.254.169.254/latest/meta-data/iam/security-credentials/').readLines()
            if( role.size() != 1 )
                throw new IllegalArgumentException("Not a valid EC2 IAM role")
            return role.get(0)
        }
        catch( IOException e ) {
            log.trace "Unable to fetch IAM credentials -- Cause: ${e.message}"
            return null
        }
    }

    /**
     * Fetch a remote URL resource text content
     *
     * @param path
     *      A valid http/https resource URL
     * @param timeout
     *      Max connection timeout in millis
     * @return
     *      The resource URL content
     */
    protected String getUrl(String path, int timeout=150) {
        final url = new URL(path)
        final con = url.openConnection()
        con.setConnectTimeout(timeout)
        con.setReadTimeout(timeout)
        return con.getInputStream().text.trim()
    }

    /**
     * Retrieve the AWS region from the EC2 instance metadata.
     * See http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-metadata.html
     *
     * @return
     *      The AWS region of the current EC2 instance eg. {@code eu-west-1} or
     *      {@code null} if it's not an EC2 instance.
     */
    protected String fetchRegion() {
        try {
            def zone = getUrl('http://169.254.169.254/latest/meta-data/placement/availability-zone')
            zone ? zone.substring(0,zone.length()-1) : null
        }
        catch (IOException e) {
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
    synchronized AmazonEC2Client getEc2Client() {

        if( ec2Client )
            return ec2Client

        def result = (accessKey && secretKey && sessionToken
                ? new AmazonEC2Client(new BasicSessionCredentials(accessKey, secretKey, sessionToken))
                : (accessKey && secretKey
                ? new AmazonEC2Client(new BasicAWSCredentials(accessKey, secretKey))
                : new AmazonEC2Client()))

        if( region )
            result.setRegion(getRegionObj(region))

        return result
    }



    /**
     * Gets or lazily creates an {@link AWSBatchClient} instance given the current
     * configuration parameter
     *
     * @return
     *      An {@link AWSBatchClient} instance
     */
    @Memoized
    AWSBatchClient getBatchClient() {
        def result = (accessKey && secretKey && sessionToken
                ? new AWSBatchClient(new BasicSessionCredentials(accessKey, secretKey, sessionToken))
                : (accessKey && secretKey
                ? new AWSBatchClient(new BasicAWSCredentials(accessKey, secretKey))
                : new AWSBatchClient()))

        if( region )
            result.setRegion(getRegionObj(region))

        return result
    }

    @Memoized
    AmazonECS getEcsClient() {

        final clientBuilder = AmazonECSClientBuilder .standard()
        if( region )
            clientBuilder.withRegion(region)

        final credentials = getCredentialsProvider0()
        if( credentials )
            clientBuilder.withCredentials(credentials)

        clientBuilder.build()
    }

    protected AWSCredentials getCredentials0() {
        if( !accessKey || !secretKey ) {
            return null
        }

        if( sessionToken )
            new BasicSessionCredentials(accessKey, secretKey, sessionToken)
        else
            new BasicAWSCredentials(accessKey, secretKey)
    }

    protected AWSStaticCredentialsProvider getCredentialsProvider0() {
        final creds = getCredentials0()
        if( !creds ) return null
        return new AWSStaticCredentialsProvider(creds)
    }


}
