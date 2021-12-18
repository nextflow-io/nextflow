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
import com.amazonaws.services.s3.model.GetObjectRequest;
import com.amazonaws.services.s3.model.S3Object;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.stream.LongStream;

/**
 * Implements a multipart downloader for S3
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class S3ParallelDownload {

    static final private Logger log = LoggerFactory.getLogger(S3ParallelDownload.class);

    private final AmazonS3 s3Client;
    private ThreadPoolExecutor executor;
    private ChunkBufferFactory bufferFactory;
    private static List<S3ParallelDownload> instances = new ArrayList<>(10);
    private final DownloadOpts opts;

    S3ParallelDownload(AmazonS3 client) {
        this(client, new DownloadOpts());
    }

    S3ParallelDownload(AmazonS3 client, DownloadOpts opts) {
        this.s3Client = client;
        this.opts = opts;
        this.executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(opts.numWorkers(), CustomThreadFactory.withName("S3-download"));
        int poolCapacity = (int) (opts.bufferMaxSize().toBytes() / opts.chunkSize());
        this.bufferFactory = new ChunkBufferFactory(opts.chunkSize(), poolCapacity);
        log.debug("Creating S3 download thread pool: {}; pool-capacity={}", opts, poolCapacity);
    }

    static public S3ParallelDownload create(AmazonS3 client, DownloadOpts opts) {
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

    static public void shutdown(boolean hard) {
        log.debug("Shutdown S3 downloader");
        for (S3ParallelDownload it : instances) {
            it.shutdown0(hard);
        }
        log.debug("Shutdown S3 downloader - done");
    }

    private Iterator<GetObjectRequest> prepareGetPartRequests(String bucketName, String key) {
        // Use range to download in parallel
        long size = s3Client.getObjectMetadata(bucketName, key).getContentLength();
        int numberOfParts = (int) Math.ceil((double) size / opts.chunkSize());
        return LongStream.range(0, numberOfParts)
                .mapToObj(index -> new AbstractMap.SimpleEntry<>(index * opts.chunkSize(),
                        (index + 1) * opts.chunkSize() > size ? size - 1 : (index + 1) * opts.chunkSize() - 1))
                .map(range -> new GetObjectRequest(bucketName, key).withRange(range.getKey(), range.getValue())).iterator();
    }

    public InputStream download(String bucketName, String key) {
        return new FutureInputStream(downloadEnumeration(prepareGetPartRequests(bucketName, key)));
    }

    private Iterator<Future<ChunkBuffer>> downloadEnumeration(Iterator<GetObjectRequest> parts) {

        // FIFO queue of parts downloading in parallel
        Queue<Future<ChunkBuffer>> futures = new LinkedList<>();

        // Add at least one part to download
        if (parts.hasNext()) {
            futures.add(executor.submit(downloadChunk(parts.next())));
        }

        // Add up to 'numWorkers' parts if the executor is idle
        int submitted = 1;
        while (parts.hasNext() && submitted < opts.numWorkers() && executor.getActiveCount() < opts.numWorkers()) {
            futures.add(executor.submit(downloadChunk(parts.next())));
            submitted++;
        }

        return new Iterator<Future<ChunkBuffer>>() {

            @Override
            public boolean hasNext() {
                return !futures.isEmpty() || parts.hasNext();
            }

            @Override
            public Future<ChunkBuffer> next() {

                // Add next part to download
                if (parts.hasNext()) {
                    futures.add(executor.submit(downloadChunk(parts.next())));
                }

                // Check if the executor is idle and submit one extra work to slowly increase the load
                if (executor.getActiveCount() < opts.numWorkers() && parts.hasNext()) {
                    futures.add(executor.submit(downloadChunk(parts.next())));
                }

                if (futures.isEmpty()) {
                    throw new NoSuchElementException();
                }

                try {
                    return futures.poll();
                } catch (Exception e) {
                    cleanUpAfterException();
                    throw new RuntimeException("Unable to complete parallel download. Individual part download failed.", e);
                }
            }

            private void cleanUpAfterException() {
                for (Future<?> future : futures) {
                    future.cancel(false);
                }
            }
        };
    }

    private Callable<ChunkBuffer> downloadChunk(final GetObjectRequest req) {
        // note: the use of the `index` determine the priority of the task in the thread pool
        return () -> {
            try (S3Object chunk = s3Client.getObject(req)) {
                final long start = req.getRange()[0];
                final long end = req.getRange()[1];
                final String path = "s3://" + req.getBucketName() + '/' + req.getKey();

                ChunkBuffer result = bufferFactory.create();
                try (InputStream stream = chunk.getObjectContent()) {
                    result.fill(stream);
                } catch (Throwable e) {
                    String msg = String.format("Failed to download chunk range=%s..%s; path=%s", start, end, path);
                    throw new IOException(msg, e);
                }
                log.trace("Downloaded chunk range={}..{}; path={}", start, end, path);
                // return it
                result.makeReadable();
                return result;
            }
        };
    }
}
