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
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.script.dsl.Description
import nextflow.exception.ProcessUnrecoverableException
import nextflow.util.Duration

/**
 * Model AWS Batch config settings
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsBatchConfig implements CloudTransferOptions, ConfigScope {

    public static final int DEFAULT_AWS_MAX_ATTEMPTS = 5

    @ConfigOption
    @Description("""
        The path where the AWS command line tool is installed in the host AMI.
    """)
    final String cliPath

    @ConfigOption
    @Description("""
        Delay between download attempts from S3 (default: `10 sec`).
    """)
    final Duration delayBetweenAttempts

    @ConfigOption
    @Description("""
        The AWS Batch [Execution Role](https://docs.aws.amazon.com/batch/latest/userguide/execution-IAM-role.html) ARN that needs to be used to execute the Batch Job. It is mandatory when using AWS Fargate.
    """)
    final String executionRole

    @ConfigOption
    @Description("""
        When `true`, add the `--force-glacier-transfer` flag to AWS CLI S3 download commands (default: `false`).
    """)
    final boolean forceGlacierTransfer

    @ConfigOption
    @Description("""
        The AWS Batch Job Role ARN that needs to be used to execute the Batch Job.
    """)
    final String jobRole

    @ConfigOption
    @Description("""
        The name of the logs group used by Batch Jobs (default: `/aws/batch/job`).
    """)
    final String logsGroup

    @ConfigOption
    @Description("""
        Max parallel upload/download transfer operations *per job* (default: `4`).
    """)
    final int maxParallelTransfers

    @ConfigOption
    @Description("""
        Max number of execution attempts of a job interrupted by a EC2 Spot reclaim event (default: `0`)
    """)
    final Integer maxSpotAttempts

    @ConfigOption
    @Description("""
        Max number of downloads attempts from S3 (default: `1`).
    """)
    final int maxTransferAttempts

    @ConfigOption
    @Description("""
        The compute platform type used by AWS Batch. Can be either `ec2` or `fargate`. Set to `fargate` to use [AWS Fargate](https://docs.aws.amazon.com/batch/latest/userguide/fargate.html).
    """)
    final String platformType

    @ConfigOption
    @Description("""
        The [retry mode](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-retries.html) used to handle rate-limiting by AWS APIs. Can be one of `standard`, `legacy`, `adaptive`, or `built-in` (default: `standard`).
    """)
    final String retryMode

    @ConfigOption
    @Description("""
        The scheduling priority for all tasks when using [fair-share scheduling](https://aws.amazon.com/blogs/hpc/introducing-fair-share-scheduling-for-aws-batch/) (default: `0`).
    """)
    final Integer schedulingPriority

    @ConfigOption
    @Description("""
        The share identifier for all tasks when using [fair-share scheduling](https://aws.amazon.com/blogs/hpc/introducing-fair-share-scheduling-for-aws-batch/).
    """)
    final String shareIdentifier

    @ConfigOption
    @Description("""
        When `true`, jobs that cannot be scheduled due to lack of resources or misconfiguration are terminated and handled as task failures (default: `false`).
    """)
    final boolean terminateUnschedulableJobs

    @ConfigOption
    @Description("""
        List of container mounts. Mounts can be specified as simple e.g. `/some/path` or canonical format e.g. `/host/path:/mount/path[:ro|rw]`.
    """)
    final List<String> volumes

    /**
     * The path for the `s5cmd` tool as an alternative to `aws s3` CLI to upload/download files
     */
    String s5cmdPath

    /**
     * Whenever it should use Fargate API
     */
    boolean fargateMode

    AwsBatchConfig(Map opts) {
        fargateMode = opts.platformType == 'fargate'
        cliPath = !fargateMode ? parseCliPath(opts.cliPath as String) : null
        s5cmdPath = fargateMode ? parses5cmdPath(opts.cliPath as String) : null
        maxParallelTransfers = opts.maxParallelTransfers as Integer ?: MAX_TRANSFER
        maxTransferAttempts = opts.maxTransferAttempts as Integer ?: defaultMaxTransferAttempts()
        delayBetweenAttempts = opts.delayBetweenAttempts as Duration ?: DEFAULT_DELAY_BETWEEN_ATTEMPTS
        maxSpotAttempts = opts.maxSpotAttempts!=null ? opts.maxSpotAttempts as Integer : null
        volumes = makeVols(opts.volumes)
        jobRole = opts.jobRole
        logsGroup = opts.logsGroup
        retryMode = opts.retryMode ?: 'standard'
        shareIdentifier = opts.shareIdentifier
        schedulingPriority = opts.schedulingPriority as Integer ?: 0
        executionRole = opts.executionRole
        terminateUnschedulableJobs = opts.terminateUnschedulableJobs as boolean
        forceGlacierTransfer = opts.forceGlacierTransfer as boolean
        if( retryMode == 'built-in' )
            retryMode = null // this force falling back on NF built-in retry mode instead of delegating to AWS CLI tool
        if( retryMode && retryMode !in AwsOptions.VALID_RETRY_MODES )
            log.warn "Unexpected value for 'aws.batch.retryMode' config setting - offending value: $retryMode - valid values: ${AwsOptions.VALID_RETRY_MODES.join(',')}"
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
