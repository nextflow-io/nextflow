/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.cloud.aws.config

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.cloud.CloudTransferOptions
import nextflow.cloud.aws.batch.AwsOptions
import nextflow.util.Duration

/**
 * Model AWS Batch config settings
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsBatchConfig implements CloudTransferOptions {

    public static final int DEFAULT_MAX_SPOT_ATTEMPTS = 5

    public static final int DEFAULT_AWS_MAX_ATTEMPTS = 5

    private int maxParallelTransfers = MAX_TRANSFER

    private int maxTransferAttempts = MAX_TRANSFER_ATTEMPTS

    private Duration delayBetweenAttempts = DEFAULT_DELAY_BETWEEN_ATTEMPTS

    private String cliPath

    private String retryMode

    private Integer maxSpotAttempts

    private Boolean debug

    private Boolean fetchInstanceType

    /**
     * The job role ARN that should be used
     */
    private String jobRole

    /**
     * The name of the logs group used by jobs
     */
    private String logsGroup

    /**
     * Volume mounts
     */
    private List<String> volumes

    /**
     * The share identifier for all tasks when using fair-share scheduling
     */
    private String shareIdentifier

    /*
     * only for testing
     */
    protected AwsBatchConfig() {}

    AwsBatchConfig(Map opts) {
        cliPath = opts.cliPath
        maxParallelTransfers = opts.maxParallelTransfers as Integer ?: MAX_TRANSFER
        maxTransferAttempts = opts.maxTransferAttempts as Integer ?: defaultMaxTransferAttempts()
        delayBetweenAttempts = opts.delayBetweenAttempts as Duration ?: DEFAULT_DELAY_BETWEEN_ATTEMPTS
        maxSpotAttempts = opts.maxSpotAttempts as Integer ?: DEFAULT_MAX_SPOT_ATTEMPTS
        volumes = makeVols(opts.volumes)
        jobRole = opts.jobRole
        logsGroup = opts.logsGroup
        fetchInstanceType = opts.fetchInstanceType
        retryMode = opts.retryMode ?: 'standard'
        shareIdentifier = opts.shareIdentifier
        if( retryMode == 'built-in' )
            retryMode = null // this force falling back on NF built-in retry mode instead of delegating to AWS CLI tool
        if( retryMode && retryMode !in AwsOptions.VALID_RETRY_MODES )
            log.warn "Unexpected value for 'aws.batch.retryMode' config setting - offending value: $retryMode - valid values: ${AwsOptions.VALID_RETRY_MODES.join(',')}"

    }

    // ====  getters =====

    String getCliPath() {
        return cliPath
    }

    int getMaxParallelTransfers() {
        return maxParallelTransfers
    }

    int getMaxTransferAttempts() {
        return maxTransferAttempts
    }

    Duration getDelayBetweenAttempts() {
        return delayBetweenAttempts
    }

    String getRetryMode() {
        return retryMode
    }

    Integer getMaxSpotAttempts() {
        return maxSpotAttempts
    }

    Boolean getDebug() {
        return debug
    }

    Boolean getFetchInstanceType() {
        return fetchInstanceType
    }

    String getJobRole() {
        return jobRole
    }

    String getLogsGroup() {
        return logsGroup
    }

    List<String> getVolumes() {
        return volumes
    }

    String getShareIdentifier() {
        return shareIdentifier
    }

    protected int defaultMaxTransferAttempts() {
        final env = SysEnv.get()
        return env.AWS_MAX_ATTEMPTS ? env.AWS_MAX_ATTEMPTS as int : DEFAULT_AWS_MAX_ATTEMPTS
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

    AwsBatchConfig addVolume(Path path) {
        assert path.scheme == 'file'
        if( volumes == null )
            volumes = new ArrayList(10)
        def location = path.toString()
        if( !volumes.contains(location) )
            volumes.add(location)
        return this
    }

}
