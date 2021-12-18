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
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.upplication.s3fs.util.ByteBufferInputStream;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.io.SequenceInputStream;
import java.nio.ByteBuffer;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Enumeration;
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
 * @author Jordi Deu-Pons <jordi@seqera.io>
 */
public class S3ParallelDownload {

    private static final int BUFFER_SIZE = 1024 * 5;
    private static final Logger log = LoggerFactory.getLogger(S3ParallelDownload.class);
    private static final List<S3ParallelDownload> instances = new ArrayList<>(10);

    private final AmazonS3 s3Client;
    private final ThreadPoolExecutor executor;
    private final DirectByteBufferPool bufferPool = new DirectByteBufferPool();
    private final DownloadOpts opts;

    S3ParallelDownload(AmazonS3 client) {
        this(client, new DownloadOpts());
    }

    S3ParallelDownload(AmazonS3 client, DownloadOpts opts) {
        this.s3Client = client;
        this.opts = opts;
        this.executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(opts.numWorkers(), new ThreadFactoryBuilder().setNameFormat("S3-download-%d").build());
        int poolCapacity = (int) (opts.bufferMaxSize().toBytes() / opts.chunkSize());
        log.debug("Creating S3 download thread pool: {}; pool-capacity={}", opts, poolCapacity);
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
        } finally {
            bufferPool.close();
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
        return new SequenceInputStream(
                getDownloadIterator(
                        prepareGetPartRequests(bucketName, key)
                )
        );
    }

    private Enumeration<InputStream> getDownloadIterator(Iterator<GetObjectRequest> parts) {

        // FIFO queue of parts downloading in parallel
        Queue<Future<InputStream>> futures = new LinkedList<>();

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

        return new Enumeration<InputStream>() {

            public boolean hasMoreElements() {
                return !futures.isEmpty() || parts.hasNext();
            }

            public InputStream nextElement() {

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
                    return futures.poll().get();
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

    private Iterator<GetObjectRequest> prepareGetPartRequests(String bucketName, String key) {
        // Use range to download in parallel
        long size = s3Client.getObjectMetadata(bucketName, key).getContentLength();
        int numberOfParts = (int) Math.ceil((double) size / opts.chunkSize());
        return LongStream.range(0, numberOfParts)
                .mapToObj(index -> new AbstractMap.SimpleEntry<>(index * opts.chunkSize(),
                        (index + 1) * opts.chunkSize() > size ? size - 1 : (index + 1) * opts.chunkSize() - 1))
                .map(range -> new GetObjectRequest(bucketName, key).withRange(range.getKey(), range.getValue())).iterator();
    }

    private Callable<InputStream> downloadChunk(GetObjectRequest req) {
        return () -> {
            try (S3Object chunk = s3Client.getObject(req)) {
                final long[] range = req.getRange();
                final String chunkRef = range == null ?
                        String.format("part=%s", req.getPartNumber()) :
                        String.format("range=%s..%s", range[0], range[1]);
                final String path = "s3://" + req.getBucketName() + '/' + req.getKey();
                final long chunkSize = chunk.getObjectMetadata().getContentLength();
                log.trace("Download chunk {}; path={}", chunkRef, path);

                try (InputStream stream = chunk.getObjectContent()) {
                    return bufferedStream(stream, chunkSize);
                } catch (OutOfMemoryError e) {
                    log.warn("Out of memory downloading chunk {}; path={}. Fallback to direct copy.", chunkRef, path);
                    return new DownloadOnDemandInputStream(req, s3Client);
                } catch (Throwable e) {
                    String msg = String.format("Failed to download chunk %s; path=%s", chunkRef, path);
                    throw new IOException(msg, e);
                }
            }
        };
    }

    private InputStream bufferedStream(InputStream stream, long size) throws IOException {

        // Create or reuse a byte buffer.
        ByteBuffer buffer = bufferPool.take((int) size);

        // Copy stream to bytebuffer
        int n;
        byte[] b = new byte[BUFFER_SIZE];
        while ((n = stream.read(b)) != -1) {
            buffer.put(b, 0, n);
        }

        // Prepare for reading
        buffer.flip();

        return new ByteBufferInputStream(buffer) {
            @Override
            public void close() {
                // Free the buffer to the pool
                bufferPool.give(buffer);
            }
        };
    }

}
