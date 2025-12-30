/*
 * Copyright 2020-2024, Seqera Labs
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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.script.dsl.Description

/**
 * Model AWS ECS Managed Instances executor config settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsEcsConfig implements ConfigScope {

    public static final String DEFAULT_LOGS_GROUP = '/aws/ecs/nextflow'
    public static final int DEFAULT_MAX_SPOT_ATTEMPTS = 5

    @ConfigOption
    @Description("""
        The name or ARN of the ECS cluster with Managed Instances capacity provider where tasks will be executed.
        This setting is required.
    """)
    final String cluster

    @ConfigOption
    @Description("""
        The ARN of the IAM role that Amazon ECS uses for task execution. This role is required for pulling
        container images and sending logs to CloudWatch.
    """)
    final String executionRole

    @ConfigOption
    @Description("""
        The ARN of the IAM role that tasks can use to make AWS API requests (e.g., for S3 access).
    """)
    final String taskRole

    @ConfigOption
    @Description("""
        List of VPC subnet IDs where ECS tasks will be launched. If not specified, subnets are
        auto-discovered from the default VPC.
    """)
    final List<String> subnets

    @ConfigOption
    @Description("""
        List of security group IDs to associate with ECS tasks. If not specified, the default
        security group is auto-discovered from the default VPC.
    """)
    final List<String> securityGroups

    @ConfigOption
    @Description("""
        The name of the CloudWatch Logs group where task logs will be sent (default: `/aws/ecs/nextflow`).
    """)
    final String logsGroup

    @ConfigOption
    @Description("""
        Maximum number of retry attempts for tasks interrupted by Spot instance reclaim (default: `5`).
    """)
    final Integer maxSpotAttempts

    @ConfigOption
    @Description("""
        Whether to assign a public IP address to tasks (default: `true`). When enabled, tasks can
        access the internet without requiring a NAT gateway.
    """)
    final Boolean assignPublicIp

    AwsEcsConfig(Map opts) {
        cluster = opts.cluster as String
        executionRole = opts.executionRole as String
        taskRole = opts.taskRole as String
        subnets = parseStringList(opts.subnets)
        securityGroups = parseStringList(opts.securityGroups)
        logsGroup = opts.logsGroup as String ?: DEFAULT_LOGS_GROUP
        maxSpotAttempts = opts.maxSpotAttempts != null ? opts.maxSpotAttempts as Integer : DEFAULT_MAX_SPOT_ATTEMPTS
        assignPublicIp = opts.assignPublicIp != null ? opts.assignPublicIp as Boolean : true
    }

    protected List<String> parseStringList(Object obj) {
        if (!obj)
            return null
        if (obj instanceof List)
            return ((List) obj).collect { it.toString() }
        if (obj instanceof CharSequence)
            return obj.toString().tokenize(',').collect { it.trim() }
        throw new IllegalArgumentException("Not a valid list value: $obj [${obj.getClass().getName()}]")
    }
}
