/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cloud.aws.batch.model

import groovy.transform.CompileStatic
import software.amazon.awssdk.services.batch.model.ContainerProperties
import software.amazon.awssdk.services.batch.model.EphemeralStorage
import software.amazon.awssdk.services.batch.model.KeyValuePair
import software.amazon.awssdk.services.batch.model.LinuxParameters
import software.amazon.awssdk.services.batch.model.LogConfiguration
import software.amazon.awssdk.services.batch.model.MountPoint
import software.amazon.awssdk.services.batch.model.NetworkConfiguration
import software.amazon.awssdk.services.batch.model.ResourceRequirement
import software.amazon.awssdk.services.batch.model.RuntimePlatform
import software.amazon.awssdk.services.batch.model.Secret
import software.amazon.awssdk.services.batch.model.Ulimit
import software.amazon.awssdk.services.batch.model.Volume

/**
 * Models the container properties used to configure an AWS Batch job.
 *
 * This is a mutable version of {@link ContainerProperties} required
 * to simplify the extension of container settings in the AWS Batch executor
 * and its sub-classes (e.g. nf-xpack).
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ContainerPropertiesModel {

    private String image

    private List<String> command

    private List<ResourceRequirement> resourceRequirements

    private String jobRoleArn

    private String executionRoleArn

    private LinuxParameters linuxParameters

    private ArrayList<KeyValuePair> environment

    private boolean privileged

    private String user

    private boolean readonlyRootFilesystem

    private ArrayList<Ulimit> ulimits

    private LogConfiguration logConfiguration

    private List<MountPoint> mountPoints

    private List<Volume> volumes

    private NetworkConfiguration networkConfiguration

    private EphemeralStorage ephemeralStorage

    private RuntimePlatform runtimePlatform

    private List<Secret> secrets

    ContainerPropertiesModel image(String value) {
        this.image = value
        return this
    }

    ContainerPropertiesModel command(String... value) {
        this.command = value as List<String>
        return this
    }

    ContainerPropertiesModel resourceRequirements(ResourceRequirement... value) {
        this.resourceRequirements = value as List<ResourceRequirement>
        return this
    }

    ContainerPropertiesModel jobRoleArn(String value) {
        this.jobRoleArn = value
        return this
    }

    ContainerPropertiesModel executionRoleArn(String value) {
        this.executionRoleArn = value
        return this
    }

    ContainerPropertiesModel user(String user) {
        this.user = user
        return this
    }

    ContainerPropertiesModel readonlyRootFilesystem(boolean value) {
        this.readonlyRootFilesystem = value
        return this
    }

    ContainerPropertiesModel environment(ArrayList<KeyValuePair> value) {
        this.environment = value
        return this
    }

    ContainerPropertiesModel linuxParameters(LinuxParameters value) {
        this.linuxParameters = value
        return this
    }

    ContainerPropertiesModel privileged(boolean value) {
        this.privileged = value
        return this
    }

    ContainerPropertiesModel ulimits(ArrayList<Ulimit> value) {
        this.ulimits = value
        return this
    }

    ContainerPropertiesModel logConfiguration(LogConfiguration value) {
        this.logConfiguration = value
        return this
    }

    ContainerPropertiesModel mountPoints(List<MountPoint> value) {
        this.mountPoints = value as List<MountPoint>
        return this
    }

    ContainerPropertiesModel volumes(List<Volume> value) {
        this.volumes = value as List<Volume>
        return this
    }

    ContainerPropertiesModel networkConfiguration(NetworkConfiguration value) {
        this.networkConfiguration = value
        return this
    }

    ContainerPropertiesModel ephemeralStorage(EphemeralStorage value) {
        this.ephemeralStorage = value
        return this
    }

    ContainerPropertiesModel runtimePlatform(RuntimePlatform value) {
        this.runtimePlatform = value
        return this
    }

    ContainerPropertiesModel secrets(List<Secret> value) {
        this.secrets = value
        return this
    }

    LinuxParameters getLinuxParameters() {
        return linuxParameters
    }

    ArrayList<KeyValuePair> getEnvironment() {
        return environment
    }

    boolean getPrivileged() {
        return privileged
    }

    String getUser() {
        return user
    }

    boolean getReadonlyRootFilesystem() {
        return readonlyRootFilesystem
    }

    ArrayList<Ulimit> getUlimits() {
        return ulimits
    }

    String getImage() {
        return image
    }

    List<String> getCommand() {
        return command
    }

    List<ResourceRequirement> getResourceRequirements() {
        return resourceRequirements
    }

    String getJobRoleArn() {
        return jobRoleArn
    }

    String getExecutionRoleArn() {
        return executionRoleArn
    }

    LogConfiguration getLogConfiguration() {
        return logConfiguration
    }

    List<MountPoint> getMountPoints() {
        return mountPoints
    }

    List<Volume> getVolumes() {
        return volumes
    }

    NetworkConfiguration getNetworkConfiguration() {
        return networkConfiguration
    }

    EphemeralStorage getEphemeralStorage() {
        return ephemeralStorage
    }

    RuntimePlatform getRuntimePlatform() {
        return runtimePlatform
    }

    ContainerProperties toBatchContainerProperties() {
        def builder = ContainerProperties.builder()

        if (image) builder.image(image)
        if (command) builder.command(command)
        if (resourceRequirements) builder.resourceRequirements(resourceRequirements)
        if (jobRoleArn) builder.jobRoleArn(jobRoleArn)
        if (executionRoleArn) builder.executionRoleArn(executionRoleArn)
        if (linuxParameters) builder.linuxParameters(linuxParameters)
        if (environment) builder.environment(environment)
        if (privileged) builder.privileged(privileged)
        if (user) builder.user(user)
        if (readonlyRootFilesystem) builder.readonlyRootFilesystem(readonlyRootFilesystem)
        if (ulimits) builder.ulimits(ulimits)
        if (logConfiguration) builder.logConfiguration(logConfiguration)
        if (mountPoints) builder.mountPoints(mountPoints)
        if (volumes) builder.volumes(volumes)
        if (networkConfiguration) builder.networkConfiguration(networkConfiguration)
        if (ephemeralStorage) builder.ephemeralStorage(ephemeralStorage)
        if (runtimePlatform) builder.runtimePlatform(runtimePlatform)
        if (secrets) builder.secrets(secrets)

        return builder.build()
    }

    @Override
    public String toString() {
        return "ContainerPropertiesModel{" +
            "image='" + image + '\'' +
            ", command=" + command +
            ", resourceRequirements=" + resourceRequirements +
            ", jobRoleArn='" + jobRoleArn + '\'' +
            ", executionRoleArn='" + executionRoleArn + '\'' +
            ", linuxParameters=" + linuxParameters +
            ", environment=" + environment +
            ", privileged=" + privileged +
            ", user='" + user + '\'' +
            ", readonlyRootFilesystem=" + readonlyRootFilesystem +
            ", ulimits=" + ulimits +
            ", logConfiguration=" + logConfiguration +
            ", mountPoints=" + mountPoints +
            ", volumes=" + volumes +
            ", networkConfiguration=" + networkConfiguration +
            ", ephemeralStorage=" + ephemeralStorage +
            ", runtimePlatform=" + runtimePlatform +
            ", secrets=" + secrets +
            '}';
    }
}
