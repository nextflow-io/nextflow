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

package nextflow.cloud.aws.ecs

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.cloud.aws.config.AwsConfig
import nextflow.cloud.aws.config.AwsEcsConfig
import nextflow.util.TestOnly

/**
 * Helper class wrapping AWS config options required for ECS Managed Instances executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
@CompileStatic
class AwsEcsOptions {

    private AwsConfig awsConfig

    @TestOnly
    protected AwsEcsOptions() {
        this.awsConfig = new AwsConfig(Collections.emptyMap())
    }

    AwsEcsOptions(AwsEcsExecutor executor) {
        awsConfig = new AwsConfig(executor.session.config.aws as Map ?: Collections.emptyMap())
    }

    AwsEcsOptions(Session session) {
        awsConfig = new AwsConfig(session.config.aws as Map ?: Collections.emptyMap())
    }

    String getRegion() {
        return awsConfig.getRegion()
    }

    String getCluster() {
        return awsConfig.ecsConfig.getCluster()
    }

    String getExecutionRole() {
        return awsConfig.ecsConfig.getExecutionRole()
    }

    String getTaskRole() {
        return awsConfig.ecsConfig.getTaskRole()
    }

    List<String> getSubnets() {
        return awsConfig.ecsConfig.getSubnets()
    }

    List<String> getSecurityGroups() {
        return awsConfig.ecsConfig.getSecurityGroups()
    }

    String getLogsGroup() {
        return awsConfig.ecsConfig.getLogsGroup()
    }

    Integer getMaxSpotAttempts() {
        return awsConfig.ecsConfig.getMaxSpotAttempts()
    }

    Boolean getAssignPublicIp() {
        return awsConfig.ecsConfig.getAssignPublicIp()
    }
}
