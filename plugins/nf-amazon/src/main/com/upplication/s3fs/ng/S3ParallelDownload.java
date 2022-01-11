/*
 * Copyright 2020-2022, Seqera Labs
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
import java.io.InterruptedIOException;
import java.net.SocketException;
import java.time.temporal.ChronoUnit;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;

import dev.failsafe.Failsafe;
import dev.failsafe.RetryPolicy;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import software.amazon.awssdk.auth.credentials.AwsBasicCredentials;
import software.amazon.awssdk.auth.credentials.AwsCredentials;
import software.amazon.awssdk.auth.credentials.AwsSessionCredentials;
import software.amazon.awssdk.auth.credentials.StaticCredentialsProvider;
import software.amazon.awssdk.core.ResponseInputStream;
import software.amazon.awssdk.regions.Region;
import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.S3ClientBuilder;
import software.amazon.awssdk.services.s3.model.GetObjectRequest;
import software.amazon.awssdk.services.s3.model.GetObjectResponse;
import software.amazon.awssdk.services.s3.model.HeadObjectRequest;

/**
 * Implements a multipart downloader for S3
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class S3ParallelDownload {

    static final private Logger log = LoggerFactory.getLogger(S3ParallelDownload.class);

    private final S3Client s3Client;
    private ExecutorService executor;
    private ChunkBufferFactory bufferFactory;
    private static List<S3ParallelDownload> instances = new ArrayList<>(10);
    private final DownloadOpts opts;

    S3ParallelDownload(S3Client client) {
        this(client, new DownloadOpts());
    }

    S3ParallelDownload(S3Client client, DownloadOpts opts) {
        if( opts.chunkSize() > opts.bufferMaxSize().toBytes() ) {
            String msg = String.format("S3 download chunk size cannot be greater than download max buffer size - offending values chunk size=%s, buffer size=%s", opts.chunkSizeMem(), opts.bufferMaxSize());
            throw new IllegalArgumentException(msg);
        }
        this.s3Client = client;
        this.opts = opts;
        this.executor = Executors.newFixedThreadPool(opts.numWorkers(), CustomThreadFactory.withName("S3-download"));
        int poolCapacity = (int)Math.ceil((float)opts.bufferMaxSize().toBytes() / opts.chunkSize());
        this.bufferFactory = new ChunkBufferFactory(opts.chunkSize(), poolCapacity);
        log.debug("Creating S3 download thread pool: {}; pool-capacity={}", opts, poolCapacity);
    }

    static public S3ParallelDownload create(Object accessKey, Object secretKey, Object sessionToken, DownloadOpts opts) {
        S3ClientBuilder builder = S3Client.builder();
        if (accessKey != null && secretKey != null) {
            AwsCredentials credentials = (sessionToken == null
                    ? AwsBasicCredentials.create(accessKey.toString(), secretKey.toString())
                    : AwsSessionCredentials.create(accessKey.toString(), secretKey.toString(), sessionToken.toString()) );

            if( System.getenv("NXF_AWS_REGION")!=null ) {
                log.debug("Creating AWS S3 client with region: " + System.getenv("NXF_AWS_REGION") );
                builder.region(Region.of(System.getenv("NXF_AWS_REGION")));
            }

            builder.credentialsProvider(StaticCredentialsProvider.create(credentials));
        }
        S3ParallelDownload result = new S3ParallelDownload(builder.build(), opts);
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

    protected List<GetObjectRequest> prepareGetPartRequests(String bucketName, String key) {
        // Use range to download in parallel
        long size = s3Client.headObject(HeadObjectRequest.builder().bucket(bucketName).key(key).build()).contentLength();
        int numberOfParts = (int) Math.ceil((double) size / opts.chunkSize());
        List<GetObjectRequest> result = new ArrayList<>(numberOfParts);
        for( int index=0; index<numberOfParts; index++ ) {
            long x = (long)index * opts.chunkSize();
            long y = (long)(index + 1) * opts.chunkSize() > size ? size - 1 : (long)(index + 1) * opts.chunkSize() - 1;
            result.add( GetObjectRequest.builder().bucket(bucketName).key(key).range(String.format("bytes=%d-%d", x, y)).build());
        }
        return result;
    }

    public InputStream download(String bucketName, String key) {
        List<GetObjectRequest> parts = prepareGetPartRequests(bucketName, key);
        Function<GetObjectRequest, ChunkBuffer> task = this::safeDownload;
        FutureIterator<GetObjectRequest, ChunkBuffer> itr = new FutureIterator<>(parts, task, executor, opts.numWorkers() * 2);
        return new FutureInputStream(itr);
    }

    private ChunkBuffer safeDownload(final GetObjectRequest req) {
        RetryPolicy<Object> retryPolicy = RetryPolicy.builder()
                .handle(SocketException.class)
                .withBackoff(50, opts.maxDelayMillis(), ChronoUnit.MILLIS)
                .withMaxAttempts(opts.maxAttempts())
                .onFailedAttempt(e -> log.error(String.format("Failed to download #%s file s3://%s/%s", req.range(), req.bucket(), req.key()), e.getLastFailure()))
                .build();

        return Failsafe.with(retryPolicy).get(() -> doDownload(req));
    }

    private ChunkBuffer doDownload(final GetObjectRequest req) throws IOException {
        try (ResponseInputStream<GetObjectResponse> stream = s3Client.getObject(req)) {
            final String path = "s3://" + req.bucket() + '/' + req.key();
            final String range = req.range();

            ChunkBuffer result = bufferFactory.create();
            try {
                result.fill(stream);
            }
            catch (Throwable e) {
                String msg = String.format("Failed to download %s path=%s", range, path);
                throw new IOException(msg, e);
            }
            log.trace("Downloaded {} path={}", range, path);
            // return it
            result.makeReadable();
            return result;
        }
        catch (InterruptedException e) {
            throw new InterruptedIOException();
        }
    }
}
