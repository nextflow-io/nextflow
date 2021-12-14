/*
 * Copyright 2020-2021, Seqera Labs
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

package com.upplication.s3fs.ng;

import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.internal.ServiceUtils;
import com.amazonaws.services.s3.model.GetObjectRequest;
import com.amazonaws.services.s3.model.S3ObjectInputStream;
import com.amazonaws.services.s3.transfer.internal.TransferManagerUtils;
import com.amazonaws.util.IOUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.SequenceInputStream;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;
import java.util.stream.LongStream;

/**
 * Implements a multipart downloader for S3
 *
 * @author Jordi Deu-Pons <jordi@seqera.io>
 */
public class S3ParallelDownload {

    private static final Logger log = LoggerFactory.getLogger(S3ParallelDownload.class);

    private final AmazonS3 s3Client;
    private final ThreadPoolExecutor executor;
    private static final List<S3ParallelDownload> instances = new ArrayList<>(10);
    private final DownloadOpts opts;

    S3ParallelDownload(AmazonS3 client) {
        this(client, new DownloadOpts());
    }

    S3ParallelDownload(AmazonS3 client, DownloadOpts opts) {
        this.s3Client = client;
        this.opts = opts;
        this.executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(opts.numWorkers());
    }

    public static S3ParallelDownload create(AmazonS3 client, DownloadOpts opts) {
        S3ParallelDownload result = new S3ParallelDownload(client, opts);
        instances.add(result);
        return result;
    }

    private void shutdown0(boolean hard) {
        if (hard)
            executor.shutdownNow();
        else
            executor.shutdown();
        try {
            executor.awaitTermination(1, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
    }

    public static void shutdown(boolean hard) {
        log.debug("Shutdown S3 downloader");
        for (S3ParallelDownload it : instances) {
            it.shutdown0(hard);
        }
        log.debug("Shutdown S3 downloader - done");
    }

    public InputStream download(String bucketName, String key) {
        return new SequenceInputStream(getParallelEnumeration(getPartRequestsIterator(bucketName, key)));
    }

    private boolean isDownloadParallel(String bucketName, String key) {
        GetObjectRequest getObjectRequest = new GetObjectRequest(bucketName, key);
        return TransferManagerUtils.isDownloadParallelizable(s3Client, getObjectRequest, ServiceUtils.getPartCount(getObjectRequest, s3Client));
    }

    private Iterator<GetObjectRequest> getPartRequestsIterator(String bucketName, String key) {
        if (isDownloadParallel(bucketName, key)) {
            // Multi-part object that can be downloaded in parallel
            int totalParts = ServiceUtils.getPartCount(new GetObjectRequest(bucketName, key), s3Client);
            return IntStream.range(1, totalParts + 1).mapToObj(p -> new GetObjectRequest(bucketName, key).withPartNumber(p)).iterator();
        } else {
            // Use range to download in parallel
            long size = s3Client.getObjectMetadata(bucketName, key).getContentLength();
            int numberOfParts = (int) Math.ceil((double) size / opts.chunkSize());
            return LongStream.range(0, numberOfParts)
                    .mapToObj(index -> new AbstractMap.SimpleEntry<>(index * opts.chunkSize(),
                            (index + 1) * opts.chunkSize() > size ? size - 1 : (index + 1) * opts.chunkSize() - 1))
                    .map(range -> new GetObjectRequest(bucketName, key).withRange(range.getKey(), range.getValue())).iterator();
        }
    }

    private Enumeration<ByteArrayInputStream> getParallelEnumeration(Iterator<GetObjectRequest> parts) {

        // FIFO queue of parts downloading in parallel
        Queue<Future<ByteArrayInputStream>> futures = new LinkedList<>();

        // Start downloading first 'numWorkers' parts
        int submitted = 0;
        while (parts.hasNext() && submitted < opts.numWorkers()) {
            GetObjectRequest nextPart = parts.next();
            futures.add(executor.submit(() -> downloadPartWithRetry(nextPart, 500, 3)));
            submitted++;
        }

        return new Enumeration<ByteArrayInputStream>() {

            @Override
            public boolean hasMoreElements() {
                return !futures.isEmpty() || parts.hasNext();
            }

            @Override
            public ByteArrayInputStream nextElement() {

                // Add next part to download
                if (parts.hasNext()) {
                    futures.add(executor.submit(() -> downloadPartWithRetry(parts.next(), 500, 3)));
                }

                if (futures.isEmpty()) {
                    throw new NoSuchElementException();
                }

                try {
                    return futures.poll().get();
                } catch (Exception e) {
                    cleanUpAfterException();
                    throw new RuntimeException("Unable to complete multipart download. Individual part download failed.", e);
                }
            }

            private void cleanUpAfterException() {
                for (Future<ByteArrayInputStream> future : futures) {
                    future.cancel(false);
                }
            }
        };
    }

    private ByteArrayInputStream downloadPart(GetObjectRequest getRequest) throws IOException {
        try (S3ObjectInputStream nextStream = s3Client.getObject(getRequest).getObjectContent()) {
            byte[] byteArray = IOUtils.toByteArray(nextStream);
            return new ByteArrayInputStream(byteArray);
        }
    }

    private ByteArrayInputStream downloadPartWithRetry(GetObjectRequest getRequest, long waitTimeMs, int maxRetries) {
        while (maxRetries > 0) {
            maxRetries--;
            try {
                return downloadPart(getRequest);
            } catch (Exception e) {
                if (maxRetries == 0) {
                    try {
                        throw e;
                    } catch (Exception ignored) { // can't happen but just in case we wrap it in
                        throw new RuntimeException(e);
                    }
                }

                log.warn("Attempt " + maxRetries + " failed", e);
                try {
                    Thread.sleep(waitTimeMs);
                } catch (InterruptedException ignored) {
                }
            }
        }

        throw new IllegalStateException();
    }

}
