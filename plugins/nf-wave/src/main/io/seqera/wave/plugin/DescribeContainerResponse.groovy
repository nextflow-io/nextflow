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

import java.time.Instant

import groovy.transform.Canonical
import groovy.transform.CompileStatic
/**
 * Model a container request record
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
@CompileStatic
class DescribeContainerResponse {

    static class User {
        Long id
        String userName
        String email
    }

    @Canonical
    static class RequestInfo {
        final User user
        final Long workspaceId
        final String containerImage
        final ContainerConfig containerConfig
        final String platform
        final String towerEndpoint
        final String fingerprint
        final Instant timestamp
        final String zoneId
        final String ipAddress
    }

    @Canonical
    static class BuildInfo {
        final String containerFile
        final String condaFile
        final String buildRepository
        final String cacheRepository
    }

    @Canonical
    static class ContainerInfo {
        String image
        String digest
    }

    final String token
    final Instant expiration
    final RequestInfo request
    final BuildInfo build
    final ContainerInfo source
    final ContainerInfo wave

}
