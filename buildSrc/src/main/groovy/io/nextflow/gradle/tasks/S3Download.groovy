package io.nextflow.gradle.tasks

import java.text.DecimalFormat

import com.amazonaws.event.ProgressEvent
import com.amazonaws.event.ProgressListener
import com.amazonaws.services.s3.transfer.Transfer
import com.amazonaws.services.s3.transfer.TransferManager
import groovy.transform.CompileStatic
import io.nextflow.gradle.tasks.AbstractS3Task
import org.gradle.api.tasks.Input
import org.gradle.api.tasks.Optional
import org.gradle.api.tasks.TaskAction

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
        TransferManager tm = new TransferManager(getS3Client())
        Transfer transfer

        // directory download
        if (keyPrefix != null) {
            logger.quiet("S3 Download recursive s3://${bucket}/${keyPrefix} → ${project.file(destDir)}/")
            transfer = (Transfer) tm.downloadDirectory(bucket, keyPrefix, project.file(destDir))
        }

        // single file download
        else {
            logger.quiet("S3 Download s3://${bucket}/${key} → ${file}")
            File f = new File(file)
            f.parentFile.mkdirs()
            transfer = (Transfer) tm.download(bucket, key, f)
        }

        S3Listener listener = new S3Listener()
        listener.transfer = transfer
        transfer.addProgressListener(listener)
        transfer.waitForCompletion()
    }

    class S3Listener implements ProgressListener {
        Transfer transfer

        DecimalFormat df = new DecimalFormat("#0.0")
        public void progressChanged(ProgressEvent e) {
            logger.info("${df.format(transfer.progress.percentTransferred)}%")
        }
    }
}
