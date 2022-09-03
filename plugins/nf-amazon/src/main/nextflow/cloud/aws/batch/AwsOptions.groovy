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

package nextflow.cloud.aws.batch

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.cloud.CloudTransferOptions
import nextflow.exception.ProcessUnrecoverableException
import nextflow.util.Duration

/**
 * Helper class wrapping AWS config options required for Batch job executions
 */
@Slf4j
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
@CompileStatic
class AwsOptions implements CloudTransferOptions {

    public static final List<String> VALID_RETRY_MODES = ['legacy','standard','adaptive']

    public static final int DEFAULT_AWS_MAX_ATTEMPTS = 5

    public static final int DEFAULT_MAX_SPOT_ATTEMPTS = 5

    private Map<String,String> env = System.getenv()

    String cliPath

    String storageClass

    String storageEncryption

    String storageKmsKeyId

    String remoteBinDir

    String region

    int maxParallelTransfers = MAX_TRANSFER

    int maxTransferAttempts = MAX_TRANSFER_ATTEMPTS

    Duration delayBetweenAttempts = DEFAULT_DELAY_BETWEEN_ATTEMPTS

    String retryMode

    int maxSpotAttempts

    Boolean debug

    volatile Boolean fetchInstanceType

    /**
     * The job role ARN that should be used
     */
    String jobRole

    /**
     * The name of the logs group used by jobs
     */
    String logsGroup

    /**
     * Volume mounts
     */
    List<String> volumes

    /**
     * The share identifier for all tasks when using fair-share scheduling
     */
    String shareIdentifier

    /**
     * @return A list of volume mounts using the docker cli convention ie. `/some/path` or `/some/path:/container/path` or `/some/path:/container/path:ro`
     */
    List<String> getVolumes() { volumes != null ? Collections.unmodifiableList(volumes) : Collections.<String>emptyList() }

    /* Only for testing purpose */
    protected AwsOptions() { }

    AwsOptions( AwsBatchExecutor executor ) {
        this(executor.session)
        this.remoteBinDir = executor.getRemoteBinDir()
    }

    AwsOptions(Session session) {
        cliPath = getCliPath0(session)
        debug = session.config.navigate('aws.client.debug') as Boolean
        storageClass = session.config.navigate('aws.client.uploadStorageClass') as String
        storageKmsKeyId = session.config.navigate('aws.client.storageKmsKeyId') as String
        storageEncryption = session.config.navigate('aws.client.storageEncryption') as String
        maxParallelTransfers = session.config.navigate('aws.batch.maxParallelTransfers', MAX_TRANSFER) as int
        maxTransferAttempts = session.config.navigate('aws.batch.maxTransferAttempts', defaultMaxTransferAttempts()) as int
        delayBetweenAttempts = session.config.navigate('aws.batch.delayBetweenAttempts', DEFAULT_DELAY_BETWEEN_ATTEMPTS) as Duration
        maxSpotAttempts = session.config.navigate('aws.batch.maxSpotAttempts', DEFAULT_MAX_SPOT_ATTEMPTS) as int
        region = session.config.navigate('aws.region') as String
        volumes = makeVols(session.config.navigate('aws.batch.volumes'))
        jobRole = session.config.navigate('aws.batch.jobRole')
        logsGroup = session.config.navigate('aws.batch.logsGroup')
        fetchInstanceType = session.config.navigate('aws.batch.fetchInstanceType')
        retryMode = session.config.navigate('aws.batch.retryMode', 'standard')
        shareIdentifier = session.config.navigate('aws.batch.shareIdentifier')
        if( retryMode == 'built-in' )
            retryMode = null // this force falling back on NF built-in retry mode instead of delegating to AWS CLI tool
        if( retryMode && retryMode !in VALID_RETRY_MODES )
            log.warn "Unexpected value for 'aws.batch.retryMode' config setting - offending value: $retryMode - valid values: ${VALID_RETRY_MODES.join(',')}"
        if( fetchInstanceType==null )
            fetchInstanceType = session.config.navigate('tower.enabled',false)
    }

    protected int defaultMaxTransferAttempts() {
        return env.AWS_MAX_ATTEMPTS ? env.AWS_MAX_ATTEMPTS as int : DEFAULT_AWS_MAX_ATTEMPTS
    }

    protected String getCliPath0(Session session) {
        def result = session.config.navigate('aws.batch.cliPath')
        if( result )
            return result

        result = session.getExecConfigProp('awsbatch','awscli',null) as String
        if( result ) {
            log.warn "Config setting `executor.awscli` has been deprecated -- Use instead `aws.batch.cliPath`"
        }
        return result
    }

    void setStorageClass(String value) {
        if( value in [null, 'STANDARD', 'STANDARD_IA', 'ONEZONE_IA', 'INTELLIGENT_TIERING', 'REDUCED_REDUNDANCY' ]) {
            this.storageClass = value
            if (value == 'REDUCED_REDUNDANCY') {
                log.warn "AWS S3 Storage Class `REDUCED_REDUNDANCY` is deprecated (and more expensive than `STANDARD`). For cost savings, look to `STANDARD_IA`, `ONEZONE_IA`, `INTELLIGENT_TIERING`."
            }
        } else {
            log.warn "Unsupported AWS storage-class: $value"
        }
    }

    void setStorageEncryption(String value) {
        if( value in [null,'AES256','aws:kms'] )
            this.storageEncryption = value
        else
            log.warn "Unsupported AWS storage-encryption: $value"
    }

    void setStorageKmsKeyId(String value) {
        this.storageKmsKeyId = value
    }

    void setCliPath(String value) {
        if( !value )
            this.cliPath = null
        else {
            if( !value.startsWith('/') ) throw new ProcessUnrecoverableException("Not a valid aws-cli tools path: $value -- it must be an absolute path")
            if( !value.endsWith('/bin/aws')) throw new ProcessUnrecoverableException("Not a valid aws-cli tools path: $value -- it must end with the `/bin/aws` suffix")
            this.cliPath = value
        }
    }

    String getAwsCli() {
        def result = getCliPath()
        if( !result ) result = 'aws'
        if( region ) result += " --region $region"
        return result
    }

    protected List<String> makeVols(obj) {
        if( !obj )
            return new ArrayList<String>(10)
        if( obj instanceof List )
            return ((List)obj).collect { normPath0(it.toString()) }
        if( obj instanceof CharSequence )
            return obj.toString().tokenize(',').collect { normPath0(it) }
        throw new IllegalArgumentException("Not a valid `aws.batch.volumes` value: $obj [${obj.getClass().getName()}]")
    }

    protected String normPath0(String it) {
        def result = it.trim()
        while( result.endsWith('/') && result.size()>1 )
            result = result.substring(0,result.size()-1)
        return result
    }

    AwsOptions addVolume(Path path) {
        assert path.scheme == 'file'
        if( volumes == null )
            volumes = new ArrayList(10)
        def location = path.toString()
        if( !volumes.contains(location) )
            volumes.add(location)
        return this
    }

}
