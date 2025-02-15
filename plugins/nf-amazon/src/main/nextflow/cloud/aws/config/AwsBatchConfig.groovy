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
import nextflow.exception.ProcessUnrecoverableException
import nextflow.util.Duration

/**
 * Model AWS Batch config settings
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsBatchConfig implements CloudTransferOptions {

    public static final int DEFAULT_MAX_SPOT_ATTEMPTS = 0

    public static final int DEFAULT_AWS_MAX_ATTEMPTS = 5

    private int maxParallelTransfers = MAX_TRANSFER

    private int maxTransferAttempts = MAX_TRANSFER_ATTEMPTS

    private Duration delayBetweenAttempts = DEFAULT_DELAY_BETWEEN_ATTEMPTS

    private String cliPath

    private String retryMode

    private Integer maxSpotAttempts

    private Boolean debug

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

    /**
     * The scheduling priority for all tasks when using fair-share scheduling (0 to 9999)
     */
    private Integer schedulingPriority

    /**
     * The container execution role
     */
    String executionRole

    /**
     * The path for the `s5cmd` tool as an alternative to `aws s3` CLI to upload/download files
     */
    String s5cmdPath

    /**
     * Whenever it should use Fargate API
     */
    boolean fargateMode

    /*
     * only for testing
     */
    protected AwsBatchConfig() {}

    AwsBatchConfig(Map opts) {
        fargateMode = opts.platformType == 'fargate'
        cliPath = !fargateMode ? parseCliPath(opts.cliPath as String) : null
        s5cmdPath = fargateMode ? parses5cmdPath(opts.cliPath as String) : null
        maxParallelTransfers = opts.maxParallelTransfers as Integer ?: MAX_TRANSFER
        maxTransferAttempts = opts.maxTransferAttempts as Integer ?: defaultMaxTransferAttempts()
        delayBetweenAttempts = opts.delayBetweenAttempts as Duration ?: DEFAULT_DELAY_BETWEEN_ATTEMPTS
        maxSpotAttempts = opts.maxSpotAttempts!=null ? opts.maxSpotAttempts as Integer : DEFAULT_MAX_SPOT_ATTEMPTS
        volumes = makeVols(opts.volumes)
        jobRole = opts.jobRole
        logsGroup = opts.logsGroup
        retryMode = opts.retryMode ?: 'standard'
        shareIdentifier = opts.shareIdentifier
        schedulingPriority = opts.schedulingPriority as Integer ?: 0
        executionRole = opts.executionRole
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

    Integer getSchedulingPriority() {
        return schedulingPriority
    }

    protected int defaultMaxTransferAttempts() {
        final env = SysEnv.get()
        return env.AWS_MAX_ATTEMPTS ? env.AWS_MAX_ATTEMPTS as int : DEFAULT_AWS_MAX_ATTEMPTS
    }

    private String parseCliPath(String value) {
        if( !value )
            return null
        if( value.tokenize('/ ').contains('s5cmd') )
            return null
        if( !value.startsWith('/') )
            throw new ProcessUnrecoverableException("Not a valid aws-cli tools path: $value -- it must be an absolute path")
        if( !value.endsWith('/bin/aws'))
            throw new ProcessUnrecoverableException("Not a valid aws-cli tools path: $value -- it must end with the `/bin/aws` suffix")
        return value
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

    protected String parses5cmdPath(String value) {
        if( !value )
            return 's5cmd'
        if( value.tokenize('/ ').contains('s5cmd') )
            return value
        return 's5cmd'
    }
}
