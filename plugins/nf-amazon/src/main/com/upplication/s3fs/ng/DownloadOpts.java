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

import java.util.Collections;
import java.util.Map;
import java.util.Properties;

import nextflow.util.Duration;
import nextflow.util.MemoryUnit;

/**
 * Model S3 download options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class DownloadOpts {

    final private boolean parallelEnabled;
    private final int queueMaxSize;
    private final int numWorkers;
    private final MemoryUnit chunkSize;
    private final MemoryUnit bufferMaxSize;
    private final int maxAttempts;
    private final Duration maxDelay;

    DownloadOpts() {
        this(new Properties(), Collections.emptyMap());
    }

    DownloadOpts(Map opts) {
        this(props(opts), Collections.emptyMap());
    }

    static private Properties props(Map opts) {
        Properties result = new Properties();
        result.putAll(opts);
        return result;
    }

    DownloadOpts(Properties props, Map<String,String> env) {
        this.parallelEnabled = props.containsKey("download_parallel")
                ? Boolean.parseBoolean(props.getProperty("download_parallel")) : (env.containsKey("NXF_S3_DOWNLOAD_PARALLEL") ? Boolean.parseBoolean(env.get("NXF_S3_DOWNLOAD_PARALLEL")) : false);

        this.queueMaxSize = props.containsKey("download_queue_max_size")
                ? Integer.parseInt(props.getProperty("download_queue_max_size")) : ( env.containsKey("NXF_S3_DOWNLOAD_QUEUE_SIZE") ? Integer.parseInt(env.get("NXF_S3_DOWNLOAD_QUEUE_SIZE")) : 10_000 );

        this.numWorkers = props.containsKey("download_num_workers")
                ? Integer.parseInt(props.getProperty("download_num_workers")) : ( env.containsKey("NXF_S3_DOWNLOAD_NUM_WORKERS") ? Integer.parseInt(env.get("NXF_S3_DOWNLOAD_NUM_WORKERS")) : 10 );

        this.chunkSize = props.containsKey("download_chunk_size")
                ? MemoryUnit.of(props.getProperty("download_chunk_size")) : ( env.containsKey("NXF_S3_DOWNLOAD_CHUNK_SIZE") ? MemoryUnit.of(env.get("NXF_S3_DOWNLOAD_CHUNK_SIZE")) : MemoryUnit.of("10 MB") );

        this.bufferMaxSize = props.containsKey("download_buffer_max_size")
                ? MemoryUnit.of(props.getProperty("download_buffer_max_size")) : ( env.containsKey("NXF_S3_DOWNLOAD_BUFFER_MAX_MEM") ? MemoryUnit.of(env.get("NXF_S3_DOWNLOAD_BUFFER_MAX_MEM")) : MemoryUnit.of("1 GB") );

        this.maxAttempts = props.containsKey("download_max_attempts")
                ? Integer.parseInt(props.getProperty("download_max_attempts")) : ( env.containsKey("NXF_S3_DOWNLOAD_MAX_ATTEMPTS") ? Integer.parseInt(env.get("NXF_S3_DOWNLOAD_MAX_ATTEMPTS")) : 5 );

        this.maxDelay = props.containsKey("download_max_delay")
                ? Duration.of(props.getProperty("download_max_delay")) : ( env.containsKey("NXF_S3_DOWNLOAD_MAX_DELAY") ? Duration.of(env.get("NXF_S3_DOWNLOAD_MAX_DELAY")) : Duration.of("90s") );

    }

    static public DownloadOpts from(Properties props) {
        return from(props, System.getenv());
    }

    static public DownloadOpts from(Properties props, Map<String,String> env) {
        return new DownloadOpts(props, env);
    }

    public boolean parallelEnabled() { return parallelEnabled; }

    @Deprecated public int queueMaxSize() { return queueMaxSize; }

    public MemoryUnit chunkSizeMem() { return chunkSize; }

    public int chunkSize() { return (int)chunkSize.toBytes(); }

    public MemoryUnit bufferMaxSize() { return bufferMaxSize; }

    public int numWorkers() { return numWorkers; }

    public long maxDelayMillis() {
        return maxDelay.getMillis();
    }

    public int maxAttempts() {
        return maxAttempts;
    }

    @Override
    public String toString() {
        return String.format("workers=%s; chunkSize=%s; queueSize=%s; max-mem=%s; maxAttempts=%s; maxDelay=%s", numWorkers, chunkSize, queueMaxSize, bufferMaxSize, maxAttempts, maxDelay);
    }
}
