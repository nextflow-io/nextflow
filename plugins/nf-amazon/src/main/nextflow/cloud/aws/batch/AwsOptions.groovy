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

package nextflow.cloud.aws.batch

import java.nio.file.Path

import com.amazonaws.services.s3.model.CannedAccessControlList
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.cloud.CloudTransferOptions
import nextflow.cloud.aws.config.AwsConfig
import nextflow.util.Duration
import nextflow.util.TestOnly
/**
 * Helper class wrapping AWS config options required for Batch job executions
 */
@Slf4j
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
@CompileStatic
class AwsOptions implements CloudTransferOptions {

    public static final List<String> VALID_RETRY_MODES = ['legacy','standard','adaptive']

    private AwsConfig awsConfig

    String remoteBinDir

    volatile Boolean fetchInstanceType

    /* Only for testing purpose */
    @TestOnly
    protected AwsOptions() {
        this.awsConfig=new AwsConfig(Collections.emptyMap())
    }

    AwsOptions( AwsBatchExecutor executor ) {
        this(executor.session)
        this.remoteBinDir = executor.getRemoteBinDir()
    }

    @Deprecated
    AwsOptions(Session session) {
        awsConfig = new AwsConfig(session.config.aws as Map ?: Collections.emptyMap())
        fetchInstanceType = session.config.navigate('aws.batch.fetchInstanceType')
        if( fetchInstanceType==null )
            fetchInstanceType = session.config.navigate('tower.enabled',false)
    }

    String getRegion() {
        return awsConfig.getRegion()
    }

    String getJobRole() {
        return awsConfig.batchConfig.getJobRole()
    }

    String getLogsGroup() {
        return awsConfig.batchConfig.getLogsGroup()
    }

    String getRetryMode() {
        return awsConfig.batchConfig.getRetryMode()
    }

    String getShareIdentifier() {
        return awsConfig.batchConfig.getShareIdentifier()
    }

    Integer getSchedulingPriority() {
        return awsConfig.batchConfig.getSchedulingPriority()
    }

    String getCliPath() {
        return awsConfig.batchConfig.getCliPath()
    }

    List<String> getVolumes() {
        final result = awsConfig.batchConfig.getVolumes()
        return result != null ? Collections.unmodifiableList(result) : Collections.<String>emptyList()
    }

    int getMaxParallelTransfers() {
        return awsConfig.batchConfig.getMaxParallelTransfers()
    }

    int getMaxTransferAttempts() {
        return awsConfig.batchConfig.getMaxTransferAttempts()
    }

    int getMaxSpotAttempts() {
        return awsConfig.batchConfig.getMaxSpotAttempts()
    }

    Duration getDelayBetweenAttempts() {
        return awsConfig.batchConfig.getDelayBetweenAttempts()
    }

    String getStorageClass() {
        return awsConfig.s3Config.getStorageClass()
    }

    String getStorageEncryption() {
        return awsConfig.s3Config.getStorageEncryption()
    }

    String getStorageKmsKeyId() {
        return awsConfig.s3Config.getStorageKmsKeyId()
    }

    CannedAccessControlList getS3Acl() {
        return awsConfig.s3Config.getS3Acl()
    }

    Boolean getDebug() {
        return awsConfig.s3Config.getDebug()
    }

    Boolean getRequesterPays() {
        return awsConfig.s3Config.getRequesterPays()
    }

    String getAwsCli() {
        def result = getCliPath()
        if( !result ) result = 'aws'
        if( region ) result += " --region $region"
        return result
    }

    AwsOptions addVolume(Path path) {
        awsConfig.batchConfig.addVolume(path)
        return this
    }

    boolean isFargateMode() {
        return awsConfig.batchConfig.fargateMode
    }

    String getS5cmdPath() {
        return awsConfig.batchConfig.s5cmdPath
    }

    String getExecutionRole() {
        return awsConfig.batchConfig.getExecutionRole()
    }

}
