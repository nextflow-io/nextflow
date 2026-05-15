/*
 * Copyright 2013-2026, Seqera Labs
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
 */
package nextflow.util.checksum

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode

/**
 * Immutable value object representing a content-addressable checksum produced
 * by a cloud storage provider (e.g. S3 CRC64NVMe, GCS CRC32C, Azure MD5).
 *
 * Used by {@link ContentChecksumResolver} as the file-content contribution to
 * the task hash when {@link nextflow.util.CacheHelper.HashMode#CLOUD_HASH} is
 * active.
 */
@CompileStatic
@EqualsAndHashCode
class Checksum {

    final String algo
    final String value

    Checksum(String algo, String value) {
        if( !algo )  throw new IllegalArgumentException("Checksum.algo must not be null or blank")
        if( !value ) throw new IllegalArgumentException("Checksum.value must not be null or blank")
        this.algo = algo
        this.value = value
    }

    String toHashContribution() { "${algo}:${value}" }
}
