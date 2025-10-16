/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.cloud.aws.nio.util;

import java.util.Properties;
import java.util.concurrent.Semaphore;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import software.amazon.awssdk.transfer.s3.S3TransferManager;
import software.amazon.awssdk.transfer.s3.model.Copy;
import software.amazon.awssdk.transfer.s3.model.CopyRequest;
import software.amazon.awssdk.transfer.s3.model.DirectoryUpload;
import software.amazon.awssdk.transfer.s3.model.DownloadFileRequest;
import software.amazon.awssdk.transfer.s3.model.FileDownload;
import software.amazon.awssdk.transfer.s3.model.FileUpload;
import software.amazon.awssdk.transfer.s3.model.UploadDirectoryRequest;
import software.amazon.awssdk.transfer.s3.model.UploadFileRequest;

import static nextflow.cloud.aws.config.AwsS3Config.*;

/**
 * Extends the S3 Transfer Manager with semaphores to limit concurrent
 * transfers based on available resources.
 *
 * Copies and uploads are limited based on the `maxConnections` setting.
 *
 * Downloads are limited based on the `maxDownloadHeapMemory` setting. The
 * CRT client allocates a buffer of 10 * part size for each transfer by default.
 *
 * @see https://github.com/aws/aws-sdk-java-v2/issues/6323
 *
 * @author Jorge Ejarque (jorge.ejarque@seqera.io)
 */
public class ExtendedS3TransferManager {

    private static final Logger log = LoggerFactory.getLogger(ExtendedS3TransferManager.class);

    private S3TransferManager transferManager;
    private Semaphore semaphore;
    private long partSize;
    private Semaphore downloadSemaphore;

    public ExtendedS3TransferManager( S3TransferManager transferManager, Properties props){
        this.transferManager = transferManager;
        setDefaultSemaphore(props);
        setDownloadSemaphore(props);
    }

    private void setDefaultSemaphore(Properties props) {
        int permits = 100;
        if( props.containsKey("max_connections")) {
            permits = Integer.parseInt(props.getProperty("max_connections"));
        }
        this.semaphore = new Semaphore(permits);
    }

    private void setDownloadSemaphore(Properties props) {
        long maxBufferSize = DEFAULT_MAX_DOWNLOAD_BUFFER_SIZE;
        if( props.containsKey("max_download_heap_memory")) {
            log.trace("AWS client config - max_download_heap_memory: {}", props.getProperty("max_download_heap_memory"));
            maxBufferSize = Long.parseLong(props.getProperty("max_download_heap_memory"));
        }

        this.partSize = DEFAULT_PART_SIZE;
        if( props.containsKey("minimum_part_size")) {
            log.trace("AWS client config - minimum_part_size: {}", props.getProperty("minimum_part_size"));
            this.partSize = Long.parseLong(props.getProperty("minimum_part_size"));
        }

        int downloadPermits = (int) Math.floor((double) maxBufferSize / partSize);
        this.downloadSemaphore = new Semaphore(downloadPermits);
    }

    public FileDownload downloadFile(DownloadFileRequest request, long size) throws InterruptedException {
        int parts = estimateParts(size);
        FileDownload fileDownload;
        downloadSemaphore.acquire(parts);
        try {
            fileDownload = transferManager.downloadFile(request);
        } catch (Throwable e) {
            // Release semaphore when runtime exception during the downloadFile submission
            downloadSemaphore.release(parts);
            throw e;
        }
        // Ensure permits are always released after completion
        fileDownload
            .completionFuture()
            .whenComplete((result, error) -> downloadSemaphore.release(parts));
        return fileDownload;
    }

    protected int estimateParts(long size) {
        if (size <= 0)
            return 1;
        int parts = (int) Math.ceil((double) size / partSize);
        return Math.min(parts, DEFAULT_INIT_BUFFER_PARTS);
    }

    public FileUpload uploadFile(UploadFileRequest request) throws InterruptedException {
        FileUpload fileUpload;
        semaphore.acquire();
        try {
            fileUpload = transferManager.uploadFile(request);
        } catch (Throwable e) {
            semaphore.release();
            throw e;
        }
        fileUpload
            .completionFuture()
            .whenComplete((result, error) -> semaphore.release());
        return fileUpload;
    }

    public DirectoryUpload uploadDirectory(UploadDirectoryRequest request) throws InterruptedException {
        DirectoryUpload directoryUpload;
        semaphore.acquire();
        try {
            directoryUpload = transferManager.uploadDirectory(request);
        } catch (Throwable e) {
            semaphore.release();
            throw e;
        }
        directoryUpload
            .completionFuture()
            .whenComplete((result, error) -> semaphore.release());
        return directoryUpload;
    }

    public Copy copy(CopyRequest request) throws InterruptedException {
        Copy copy;
        semaphore.acquire();
        try {
            copy = transferManager.copy(request);
        } catch (Throwable e) {
            semaphore.release();
            throw e;
        }
        copy
            .completionFuture()
            .whenComplete((result, error) -> semaphore.release());
        return copy;
    }

}
