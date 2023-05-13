/*
 * Copyright 2013-2023, Seqera Labs
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
     * The ID of the executed container
     */
    String workflowId

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
    String condaFile

    /**
     * Spack recipe file used to build the container
     */
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

}
