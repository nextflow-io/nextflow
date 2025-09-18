/*
 * Copyright 2020-2025, Seqera Labs
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

import java.util.concurrent.ExecutionException;

import software.amazon.awssdk.transfer.s3.S3TransferManager;
import software.amazon.awssdk.transfer.s3.model.*;

/**
 * Synchronous wrapper class for the S3 Transfer Manager.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class S3TransferManagerSync {

    private S3TransferManager transferManager;

    public S3TransferManagerSync(S3TransferManager transferManager) {
        this.transferManager = transferManager;
    }

    public CompletedFileDownload downloadFile(DownloadFileRequest request) throws InterruptedException, ExecutionException {
        return transferManager.downloadFile(request).completionFuture().get();
    }

    public CompletedDirectoryDownload downloadDirectory(DownloadDirectoryRequest request) throws InterruptedException, ExecutionException {
        return transferManager.downloadDirectory(request).completionFuture().get();
    }

    public CompletedFileUpload uploadFile(UploadFileRequest request) throws InterruptedException, ExecutionException {
        return transferManager.uploadFile(request).completionFuture().get();
    }

    public CompletedDirectoryUpload uploadDirectory(UploadDirectoryRequest request) throws InterruptedException, ExecutionException {
        return transferManager.uploadDirectory(request).completionFuture().get();
    }

    public CompletedCopy copy(CopyRequest request) throws InterruptedException, ExecutionException {
        return transferManager.copy(request).completionFuture().get();
    }

}
