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
