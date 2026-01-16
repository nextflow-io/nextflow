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

package io.seqera.config

import groovy.transform.CompileStatic
import nextflow.util.Duration

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class SeqeraConfig {

    private RetryOpts retryOpts

    private String endpoint

    private String region

    private String keyPairName

    private Duration batchFlushInterval

    SeqeraConfig(Map opts) {
        this.retryOpts = new RetryOpts(opts.retryPolicy as Map ?: Map.of())
        this.endpoint = opts.endpoint as String
        if (!endpoint)
            throw new IllegalArgumentException("Missing Seqera endpoint - make sure to specify 'seqera.endpoint' settings")

        this.region = opts.region as String
        if (!region)
            region = "eu-central-1"

        this.keyPairName = opts.keyPairName as String

        // Batch submission configuration
        this.batchFlushInterval = opts.batchFlushInterval
            ? Duration.of(opts.batchFlushInterval as String)
            : Duration.of('1 sec')
    }

    RetryOpts retryOpts() {
        this.retryOpts
    }

    String getEndpoint() {
        return endpoint
    }

    String getRegion() {
        return region
    }

    String getKeyPairName() {
        return keyPairName
    }

    Duration getBatchFlushInterval() {
        return batchFlushInterval
    }
}
