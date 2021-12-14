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

import java.io.IOException;
import java.io.InputStream;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.LongStream;

import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.model.GetObjectRequest;
import com.amazonaws.services.s3.model.S3Object;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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
    private static AtomicInteger chunksCount = new AtomicInteger();
    private static AtomicInteger filesCount = new AtomicInteger();
    private final DownloadOpts opts;

    S3ParallelDownload(AmazonS3 client) {
        this(client, new DownloadOpts());
    }

    S3ParallelDownload(AmazonS3 client, DownloadOpts opts) {
        this.s3Client = client;
        this.opts = opts;
        this.executor = PriorityThreadPool.create("S3-downloader", opts.numWorkers(), opts.queueMaxSize());
        int poolCapacity = (int)(opts.bufferMaxSize().toBytes() / opts.chunkSize());
        this.bufferFactory = new ChunkBufferFactory(opts.chunkSize(), poolCapacity);
        log.debug("Creating S3 download thread pool: {}; pool-capacity={}", opts, poolCapacity);
    }

    static public S3ParallelDownload create(AmazonS3 client, DownloadOpts opts) {
        S3ParallelDownload result = new S3ParallelDownload(client, opts);
        instances.add(result);
        return result;
    }

    private void shutdown0(boolean hard) {
        if( hard )
            executor.shutdownNow();
        else
            executor.shutdown();
        try {
            executor.awaitTermination(1, TimeUnit.HOURS);
        }
        catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
    }

    static public void shutdown(boolean hard) {
        log.debug("Shutdown S3 downloader");
        for( S3ParallelDownload it : instances )  {
            it.shutdown0(hard);
        }
        log.debug("Shutdown S3 downloader - done");
    }
    
    protected List<GetObjectRequest> prepareGetPartRequests(String bucketName, String key, long size) {
        int numberOfParts = (int)Math.ceil( (double)size / opts.chunkSize() );
        return LongStream.range(0, numberOfParts)
                .mapToObj(index -> new AbstractMap.SimpleEntry<>(index * opts.chunkSize(),
                        (index + 1) * opts.chunkSize() > size ? size - 1 : (index + 1) * opts.chunkSize() - 1))
                .map(range -> new GetObjectRequest(bucketName, key).withRange(range.getKey(), range.getValue()))
                .collect(Collectors.toList());
    }

    public InputStream download( String bucketName, String key ) {
        long length = s3Client.getObjectMetadata(bucketName, key) .getContentLength();
        ChunkedInputStream result = new ChunkedInputStream(length);

        List<GetObjectRequest> chunks = prepareGetPartRequests(bucketName, key, length);
        if( length>0 && chunks.size()==0 )
            throw new IllegalStateException(String.format("Object length is greater %s but got zero chunks to download", length));

        int chunkIndex=0;
        int fileIndex = filesCount.getAndIncrement();
        for( GetObjectRequest it : chunks ) {
            executor.submit(downloadChunk(result, fileIndex, chunkIndex++, it));
            // give the ability to schedule chunks of other file transfers 
            try { Thread.sleep(10); }
            catch (InterruptedException e) { Thread.currentThread().interrupt(); }
        }

        return result;
    }

    static int seqOrder(int fileIndex, int chunkIndex) {
        return (1 << 16) * fileIndex + chunkIndex;
    }

    protected int getPriorityIndex(final int fileIndex, final int chunkIndex) {
        if( opts.strategy() == DownloadOpts.Strategy.interleaved ) {
            // give an incremental index to each across all files
            // this should behave as first arrived first served irrespective the file
            return chunksCount.incrementAndGet();
        }
        else if( opts.strategy() == DownloadOpts.Strategy.sequential ) {
            return seqOrder(fileIndex, chunkIndex);
        }
        throw new IllegalArgumentException("Unexpected download strategy=" + opts.strategy());
    }

    private Runnable downloadChunk(ChunkedInputStream chunkedStream, final int fileIndex, final int chunkIndex, final GetObjectRequest req) {
        // note: the use of the `index` determine the priority of the task in the thread pool
        return new PriorityThreadPool.PriorityRunnable(getPriorityIndex(fileIndex,chunkIndex)) {
            @Override
            public void run()  {
                try ( S3Object chunk = s3Client.getObject(req) ) {
                    final long start = req.getRange()[0];
                    final long end = req.getRange()[1];
                    final String path = "s3://" + req.getBucketName() + '/' + req.getKey();

                    ChunkBuffer result = bufferFactory.create(chunkIndex);
                    try ( InputStream stream = chunk.getObjectContent() ) {
                        result.fill(stream);
                    }
                    catch (Throwable e) {
                        String msg = String.format("Failed to download chunk index=%s; range=%s..%s; path=%s", chunkIndex, start, end, path);
                        throw new IOException(msg, e);
                    }
                    log.trace("Downloaded chunk index={}; range={}..{}; path={}", chunkIndex, start, end, path);
                    // return it
                    result.makeReadable();
                    chunkedStream.add(result);
                }
                catch (IOException e) {
                    chunkedStream.throwError(e);
                }
                catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                }
            }
        };
    }

}
