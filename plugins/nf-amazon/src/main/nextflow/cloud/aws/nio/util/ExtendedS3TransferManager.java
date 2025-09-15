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

import nextflow.cloud.aws.nio.S3Client;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import software.amazon.awssdk.transfer.s3.S3TransferManager;
import software.amazon.awssdk.transfer.s3.model.*;

import java.util.Properties;
import java.util.concurrent.Semaphore;

/**
 * Extends the S3 Transfer manager functionality limiting the number of concurrent downloads according to a target heap memory consumption.
 * Implemented to fix https://github.com/aws/aws-sdk-java-v2/issues/6323 where several concurrent downloads with CRT client can consume large heap memory space.
 * This is due to the inital download buffer allocated for each transfer that by default it is 10 * part size.
 *
 * @author Jorge Ejarque (jorge.ejarque@seqera.io)
 */
public class ExtendedS3TransferManager {

    private static final Logger log = LoggerFactory.getLogger(S3Client.class);
    private S3TransferManager transferManager;

    // According to CRT Async client docs https://sdk.amazonaws.com/java/api/latest/software/amazon/awssdk/services/s3/S3CrtAsyncClientBuilder.html
    private static final long DEFAULT_PART_SIZE = 8 * 1024 * 1024; // 8 MB
    private static final int DEFAULT_INIT_BUFFER_PARTS = 10;
    // Maximum heap buffer size
    private static final long DEFAULT_MAX_DOWNLOAD_BUFFER_SIZE = 400 * 1024 * 1024; // 8 MB
    private Semaphore concurrentDownloadSemaphore;
    private long partSize;
    private int maxPartsTotal;


    public ExtendedS3TransferManager( S3TransferManager transferManager, Properties props){
        this.transferManager = transferManager;
        setDownloadBufferProperties(props);
        this.concurrentDownloadSemaphore = new Semaphore(maxPartsTotal);
    }

    public FileDownload downloadFile(DownloadFileRequest downloadFileRequest, long size) throws InterruptedException {
        final int parts = estimateParts(size);

        concurrentDownloadSemaphore.acquire(parts);

        FileDownload downloadFile = transferManager.downloadFile(downloadFileRequest);

        // Ensure permits are always released after completion
        downloadFile.completionFuture().whenComplete((result, error) -> {
            concurrentDownloadSemaphore.release(parts); // release no matter what
        });
        return downloadFile;
    }

    private void setDownloadBufferProperties(Properties props) {
        long maxBuffer = DEFAULT_MAX_DOWNLOAD_BUFFER_SIZE;
        if( props.containsKey("max_download_buffer")) {
            log.trace("AWS client config - max_download_buffer: {}", props.getProperty("max_download_buffer"));
            maxBuffer = Long.parseLong(props.getProperty("max_download_buffer"));
        }
        this.partSize = DEFAULT_PART_SIZE;
        if( props.containsKey("minimum_part_size")) {
            log.trace("AWS client config - minimum_part_size: {}", props.getProperty("minimum_part_size"));
            this.partSize = Long.parseLong(props.getProperty("minimum_part_size"));
        }
        this.maxPartsTotal = (int) Math.floor((double) maxBuffer / partSize);
    }

    protected int estimateParts(long size) {
        if (size <= 0)
            return 1;
        int parts = (int) Math.ceil((double) size / partSize);
        return Math.min(parts, DEFAULT_INIT_BUFFER_PARTS);
    }

    public FileUpload uploadFile(UploadFileRequest request) {
        return transferManager.uploadFile(request);
    }

    public DirectoryUpload uploadDirectory(UploadDirectoryRequest request)  {
        return transferManager.uploadDirectory(request);
    }

    public Copy copy(CopyRequest request)  {
        return transferManager.copy(request);
    }

    public long getPartSize() {
        return partSize;
    }

    public int getMaxPartsTotal() {
        return maxPartsTotal;
    }


}
