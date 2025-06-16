package software.amazon.nio.spi.s3

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.aws.util.AwsHelper
import software.amazon.awssdk.core.client.config.ClientOverrideConfiguration
import software.amazon.awssdk.services.s3.model.*

import java.nio.file.Path

@Slf4j
@CompileStatic
class NextflowS3PathOpenOptions extends S3OpenOption {
    private List<Tag> tags
    private String storageClass
    private String contentType


    protected NextflowS3PathOpenOptions(List<Tag> tags, String storageClass , String contentType) {
        this.tags = tags
        this.storageClass = storageClass
        this.contentType = contentType
    }

    @Override
    S3OpenOption copy() {
        return new NextflowS3PathOpenOptions(this.tags, this.storageClass, this.contentType)
    }

    @Override
    protected void apply(GetObjectRequest.Builder getObjectRequest) {
    }

    @Override
    protected void apply(PutObjectRequest.Builder putObjectRequest, Path file) {
        if( tags )
            putObjectRequest.tagging(Tagging.builder().tagSet(tags).build())
        if( storageClass )
            putObjectRequest.storageClass(storageClass)
        if( contentType )
            putObjectRequest.contentType(contentType)

    }
}
