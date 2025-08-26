package io.nextflow.gradle.tasks

import org.gradle.api.GradleException

import java.text.DecimalFormat
import groovy.transform.CompileStatic
import io.nextflow.gradle.tasks.AbstractS3Task
import org.gradle.api.tasks.Input
import org.gradle.api.tasks.Optional
import org.gradle.api.tasks.TaskAction
import software.amazon.awssdk.transfer.s3.S3TransferManager
import software.amazon.awssdk.transfer.s3.model.DownloadDirectoryRequest
import software.amazon.awssdk.transfer.s3.model.DownloadFileRequest
import software.amazon.awssdk.transfer.s3.progress.TransferListener

import java.util.concurrent.ExecutionException

/**
 * Download files from a S3 bucket
 *
 * Based on https://github.com/mgk/s3-plugin/blob/master/src/main/groovy/com/github/mgk/gradle/S3Plugin.groovy
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class S3Download extends AbstractS3Task {

    @Input
    String bucket

    @Input
    String key

    @Input
    String file

    @Optional
    @Input
    String keyPrefix

    @Optional
    @Input
    String destDir

    @TaskAction
    def task() {
        def transferManager = S3TransferManager.builder()
            .s3Client(getS3AsyncClient()) // ensure this returns an `S3AsyncClient`
            .build()

        def listener = new S3Listener()

        // directory download
        if (keyPrefix != null) {
            logger.quiet("S3 Download recursive s3://${bucket}/${keyPrefix} → ${project.file(destDir)}/")
            def request = DownloadDirectoryRequest.builder()
                .bucket(bucket)
                .listObjectsV2RequestTransformer(builder -> builder.prefix(keyPrefix))
                .destination(project.file(destDir).toPath())
                .downloadFileRequestTransformer(builder-> builder.addTransferListener(listener))
                .build()

            def download = transferManager.downloadDirectory(request)
            try{
			    def completed = download.completionFuture().get()
                if (!completed.failedTransfers().isEmpty()){
                    throw new GradleException("Some transfers in S3 download directory: s3://"+ bucket +"/"+ keyPrefix +" has failed - Transfers: " +  completed.failedTransfers() );
			    }
		    } catch (InterruptedException e){
			    logger.debug("S3 download directory: s3://{}/{} interrupted", bucket, keyPrefix);
			    Thread.currentThread().interrupt();
		    } catch ( ExecutionException e) {
                String msg = String.format("Exception thrown downloading S3 object s3://{}/{}", bucket, keyPrefix);
                throw new GradleException(msg, e.getCause());
            }

        // single file download
        } else {
            logger.quiet("S3 Download s3://${bucket}/${key} → ${file}")
            File f = new File(file)
            f.parentFile.mkdirs()
            def request = DownloadFileRequest.builder()
                .getObjectRequest { it.bucket(bucket).key(key) }
                .destination(f.toPath())
                .addTransferListener(listener)
                .build()

            def download = transferManager.downloadFile(request)
            try{
                download.completionFuture().get()
            } catch (InterruptedException e){
			    logger.debug("S3 download file: s3://{}/{} interrupted", bucket, key);
			    Thread.currentThread().interrupt();
		    } catch ( ExecutionException e) {
                String msg = String.format("Exception thrown downloading S3 object s3://{}/{}", bucket, key);
                throw new GradleException(msg, e.getCause());
            }
        }

        transferManager.close()
    }

    class S3Listener implements TransferListener {
        DecimalFormat df = new DecimalFormat("#0.0")

        @Override
        void bytesTransferred(Context.BytesTransferred context) {
            def request = context.request() as DownloadFileRequest
            if (!context.progressSnapshot().ratioTransferred().empty)
                logger.info("${request.destination()}: ${df.format(context.progressSnapshot().ratioTransferred().asDouble)}%")
        }
    }
}
