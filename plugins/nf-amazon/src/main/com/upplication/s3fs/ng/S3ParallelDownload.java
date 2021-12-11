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
import nextflow.util.MemoryUnit;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Implements a multipart downloader for S3
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class S3ParallelDownload {

    static final private Logger log = LoggerFactory.getLogger(S3ParallelDownload.class);

    private AmazonS3 s3Client;
    private int chunkSize = chunkSize();
    private int workers = numWorkers();
    private MemoryUnit bufferMaxMem = bufferMaxSize();
    private int queueSize = 10_000;
    private ThreadPoolExecutor executor;
    private ChunkBufferFactory bufferFactory;
    private static List<S3ParallelDownload> instances = new ArrayList<>(10);
    private static AtomicInteger chunksCount = new AtomicInteger();

    private static int chunkSize() {
        if( System.getenv("NXF_S3_DOWNLOAD_CHUNK_SIZE")!=null )
            return Integer.parseInt(System.getenv("NXF_S3_DOWNLOAD_CHUNK_SIZE"));
        else
            return 10 * 1024 * 1024;
    }

    private static int numWorkers() {
        if( System.getenv("NXF_S3_DOWNLOAD_NUM_WORKERS")!=null )
            return Integer.parseInt(System.getenv("NXF_S3_DOWNLOAD_NUM_WORKERS"));
        else
            return 10;
    }

    private static MemoryUnit bufferMaxSize() {
        if( System.getenv("NXF_S3_DOWNLOAD_BUFFER_MAX_MEM")!=null )
            return MemoryUnit.of(System.getenv("NXF_S3_DOWNLOAD_BUFFER_MAX_MEM"));
        else
            return MemoryUnit.of("1 GB");
    }

    private static int queueSize() {
        if( System.getenv("NXF_S3_DOWNLOAD_QUEUE_SIZE")!=null )
            return Integer.parseInt(System.getenv("NXF_S3_DOWNLOAD_QUEUE_SIZE"));
        else
            return 10_000;
    }

    private S3ParallelDownload(AmazonS3 client) {
        this.s3Client = client;
        this.executor = PriorityThreadPool.create("S3-downloader", workers, queueSize);
        int poolCapacity = (int)(bufferMaxMem.toBytes() / chunkSize);
        this.bufferFactory = new ChunkBufferFactory(chunkSize, poolCapacity);
        log.debug("Creating S3 download thread pool: workers={}; chunkSize={}; queueSize={}; max-mem={}; buffers={}", workers, chunkSize, queueSize, bufferMaxMem, poolCapacity);
    }

    static public S3ParallelDownload create(AmazonS3 client) {
        S3ParallelDownload result = new S3ParallelDownload(client);
        instances.add(result);
        return result;
    }

    private void shutdown0() {
        executor.shutdown();
        try {
            executor.awaitTermination(1, TimeUnit.HOURS);
        }
        catch (InterruptedException e) {
            log.debug("Thread pool await termination was interrupted", e);
        }
    }

    static public void shutdown() {
        log.debug("Shutdown S3 downloader");
        for( S3ParallelDownload it : instances )  {
            it.shutdown0();
        }
        log.debug("Shutdown S3 downloader - done");
    }
    
    S3ParallelDownload withChunkSize(int chunkSize) {
        if( chunkSize<= 0 )
            throw new IllegalArgumentException(String.format("Download chunkSize cannot be less or equals to zero: %s", chunkSize));
        this.chunkSize = chunkSize;
        return this;
    }

    protected List<GetObjectRequest> prepareGetPartRequests(String bucketName, String key, long size) {
        int numberOfParts = (int)Math.ceil( (double)size / chunkSize );
        return LongStream.range(0, numberOfParts)
                .mapToObj(index -> new AbstractMap.SimpleEntry<>(index * chunkSize,
                        (index + 1) * chunkSize > size ? size - 1 : (index + 1) * chunkSize - 1))
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
        for( GetObjectRequest it : chunks ) {
            executor.submit(downloadChunk(result, chunkIndex++, it));
        }

        return result;
    }


    private Runnable downloadChunk(ChunkedInputStream chunkedStream, final int chunkIndex, final GetObjectRequest req) {
        // note: the use of the `index` determine the priority of the task in the thread pool
        return new PriorityThreadPool.PriorityRunnable(chunksCount.incrementAndGet()) {
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
                        log.debug("Failed to download chunk index={}; range={}..{}; path={}", chunkIndex, start, end, path);
                        throw e;
                    }
                    log.trace("Downloaded chunk index={}; range={}..{}; path={}", chunkIndex, start, end, path);
                    // return it
                    chunkedStream.add(result);
                }
                catch (IOException | InterruptedException e) {
                    throw new RuntimeException(e);
                }
            }
        };
    }

}
