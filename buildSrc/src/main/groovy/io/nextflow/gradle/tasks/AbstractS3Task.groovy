package io.nextflow.gradle.tasks

import org.gradle.api.DefaultTask
import org.gradle.api.tasks.Internal
import software.amazon.awssdk.auth.credentials.DefaultCredentialsProvider
import software.amazon.awssdk.services.s3.S3AsyncClient
import software.amazon.awssdk.services.s3.S3Client
import software.amazon.awssdk.regions.Region

/**
 * S3 common operations
 *
 * Based on https://github.com/mgk/s3-plugin/blob/master/src/main/groovy/com/github/mgk/gradle/S3Plugin.groovy
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
abstract class AbstractS3Task extends DefaultTask {

    @Internal
    S3Client getS3Client() {
        def credBuilder = DefaultCredentialsProvider.builder()
        if (project.s3.profile) {
            logger.quiet("Using AWS credentials profile: ${project.s3.profile}")
            credBuilder.profileName(project.s3.profile)
        }
        def clientBuilder = S3Client.builder().credentialsProvider(credBuilder.build())
        String region = project.s3.region
        if (region) {
            clientBuilder.region(Region.of(region))
        }
        return clientBuilder.build()
    }

    @Internal
    S3AsyncClient getS3AsyncClient() {
        def credBuilder = DefaultCredentialsProvider.builder()
        if (project.s3.profile) {
            logger.quiet("Using AWS credentials profile: ${project.s3.profile}")
            credBuilder.profileName(project.s3.profile)
        }
        def clientBuilder = S3AsyncClient.crtBuilder().credentialsProvider(credBuilder.build())
        String region = project.s3.region
        if (region) {
            clientBuilder.region(Region.of(region))
        }
        return clientBuilder.build()
    }
}
