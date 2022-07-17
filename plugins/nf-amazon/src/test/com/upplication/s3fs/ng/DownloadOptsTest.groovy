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

package com.upplication.s3fs.ng

import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DownloadOptsTest extends Specification {

    def 'should get default options' () {
        given:
        def props = new Properties()

        when:
        def opts = DownloadOpts.from(props)
        then:
        opts.numWorkers() == 10
        opts.queueMaxSize() == 10_000
        opts.bufferMaxSize() == MemoryUnit.of('1 GB')
        opts.chunkSize() == 10 * 1024 * 1024
        !opts.parallelEnabled()
        opts.maxDelayMillis() == Duration.of('90s').toMillis()
        opts.maxAttempts() == 5
    }

    def 'should set options with properties' () {
        given:
        def CONFIG = '''
        download_parallel = false
        download_queue_max_size = 11
        download_buffer_max_size = 222MB
        download_num_workers = 33
        download_chunk_size = 44
        download_max_attempts = 99
        download_max_delay = 99s
        '''
        def props = new Properties()
        props.load(new StringReader(CONFIG))

        when:
        def opts = DownloadOpts.from(props)
        then:
        opts.numWorkers() == 33
        opts.queueMaxSize() == 11
        opts.bufferMaxSize() == MemoryUnit.of('222 MB')
        opts.chunkSize() == 44
        !opts.parallelEnabled()
        opts.maxAttempts() == 99
        opts.maxDelayMillis() == Duration.of('99s').toMillis()
    }


    def 'should set options with env' () {
        given:
        def ENV = [
                NXF_S3_DOWNLOAD_PARALLEL: 'false',
                NXF_S3_DOWNLOAD_QUEUE_SIZE: '11',
                NXF_S3_DOWNLOAD_NUM_WORKERS: '22',
                NXF_S3_DOWNLOAD_CHUNK_SIZE: '33',
                NXF_S3_DOWNLOAD_BUFFER_MAX_MEM: '44 G',
                NXF_S3_DOWNLOAD_MAX_ATTEMPTS: '88',
                NXF_S3_DOWNLOAD_MAX_DELAY: '88s'
        ]

        when:
        def opts = DownloadOpts.from(new Properties(), ENV)
        then:
        !opts.parallelEnabled()
        opts.queueMaxSize() == 11
        opts.numWorkers() == 22
        opts.chunkSize() == 33
        opts.bufferMaxSize() == MemoryUnit.of('44 GB')
        opts.maxAttempts() == 88
        opts.maxDelayMillis() == Duration.of('88s').toMillis()
    }

}
