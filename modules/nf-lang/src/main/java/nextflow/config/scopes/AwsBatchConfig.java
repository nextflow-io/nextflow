/*
 * Copyright 2024-2025, Seqera Labs
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
package nextflow.config.scopes;

import java.util.List;

import nextflow.config.schema.ConfigOption;
import nextflow.config.schema.ConfigScope;
import nextflow.script.dsl.Description;
import nextflow.script.types.Duration;

public class AwsBatchConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The path where the AWS command line tool is installed in the host AMI.
    """)
    public String cliPath;

    @ConfigOption
    @Description("""
        Delay between download attempts from S3 (default: `10 sec`).
    """)
    public Duration delayBetweenAttempts;

    @ConfigOption
    @Description("""
        The AWS Batch Execution Role ARN that needs to be used to execute the Batch Job.

        [Read more](https://docs.aws.amazon.com/batch/latest/userguide/execution-IAM-role.html)
    """)
    public String executionRole;

    @ConfigOption
    @Description("""
        The AWS Batch Job Role ARN that needs to be used to execute the Batch Job.
    """)
    public String jobRole;

    @ConfigOption
    @Description("""
        The name of the logs group used by Batch Jobs (default: `/aws/batch`).
    """)
    public String logsGroup;

    @ConfigOption
    @Description("""
        Max parallel upload/download transfer operations *per job* (default: `4`).
    """)
    public int maxParallelTransfers;

    @ConfigOption
    @Description("""
        Max number of execution attempts of a job interrupted by a EC2 spot reclaim event (default: `5`)
    """)
    public int maxSpotAttempts;

    @ConfigOption
    @Description("""
        Max number of downloads attempts from S3 (default: `1`).
    """)
    public int maxTransferAttempts;

    @ConfigOption
    @Description("""
        The compute platform type used by AWS Batch. Can be either `ec2` or `fargate`.
    """)
    public String platformType;

    @ConfigOption
    @Description("""
        The retry mode used to accommodate rate-limiting on AWS services. Can be one of `standard`, `legacy`, `adaptive`, or `built-in` (default: `standard`).

        [Read more](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-retries.html)
    """)
    public String retryMode;

    @ConfigOption
    @Description("""
        The scheduling priority for all tasks when using fair-share scheduling for AWS Batch (default: `0`).

        [Read more](https://aws.amazon.com/blogs/hpc/introducing-fair-share-scheduling-for-aws-batch/)
    """)
    public int schedulingPriority;

    @ConfigOption
    @Description("""
        The share identifier for all tasks when using fair-share scheduling for AWS Batch.

        [Read more](https://aws.amazon.com/blogs/hpc/introducing-fair-share-scheduling-for-aws-batch/)
    """)
    public String shareIdentifier;

    @ConfigOption
    @Description("""
        One or more container mounts. Mounts can be specified as simple e.g. `/some/path` or canonical format e.g. `/host/path:/mount/path[:ro|rw]`.
    """)
    public List<String> volumes;

}
