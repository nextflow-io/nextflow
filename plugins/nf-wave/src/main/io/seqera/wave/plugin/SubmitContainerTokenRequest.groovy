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
 *
 */

package io.seqera.wave.plugin


import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import io.seqera.wave.api.ImageNameStrategy
import io.seqera.wave.api.PackagesSpec
import io.seqera.wave.api.ScanLevel
import io.seqera.wave.api.ScanMode

/**
 * Model a request for an augmented container
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@EqualsAndHashCode
@ToString(includeNames = true, includePackage = false)
@CompileStatic
class SubmitContainerTokenRequest {

    /**
     * Tower access token required to enable the service
     */
    String towerAccessToken

    /**
     * Tower refresh token
     */
    String towerRefreshToken

    /**
     * Tower workspace id
     */
    Long towerWorkspaceId

    /**
     * Tower endpoint
     */
    String towerEndpoint

    /**
     * Container image to be pulled
     */
    String containerImage

    /**
     * Container build file i.g. Dockerfile of the container to be build
     */
    String containerFile

    /**
     * List of layers to be added in the pulled image
     */
    ContainerConfig containerConfig

    /**
     * Conda recipe file used to build the container
     */
    @Deprecated
    String condaFile

    /**
     * Spack recipe file used to build the container
     */
    @Deprecated
    String spackFile

    /**
     * The request container platform
     */
    String containerPlatform

    /**
     * The target repository where the built container needs to be stored
     */
    String buildRepository

    /**
     * The container repository to cache build layers
     */
    String cacheRepository

    /**
     * Request
     */
    String timestamp

    /**
     * Request unique fingerprint
     */
    String fingerprint

    /**
     * Enable freeze container mode
     */
    boolean freeze

    /**
     * Specify the format of the container file
     */
    String format

    /**
     * When {@code true} build requests are carried out in dry-run mode.
     */
    Boolean dryRun

    /**
     * Id of compute workflow environment in tower
     */
    String workflowId

    /**
     * One or more container should be included in upstream container request
     */
    List<String> containerIncludes

    /**
     * Defines the packages to be included in this container request
     */
    PackagesSpec packages

    /**
     * The strategy applied to name a container build by wave when using
     * the freeze option.
     */
    ImageNameStrategy nameStrategy;

    /**
     * Whenever use container "mirror" mode
     */
    boolean mirror;

    /**
     * The request security scan mode
     */
    ScanMode scanMode;

    /**
     * Define the allows security vulnerabilities in the container request.
     * Empty or null means no vulnerabilities are allowed.
     */
    List<ScanLevel> scanLevels

}
