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
package nextflow.file

import java.nio.file.Path

import nextflow.util.checksum.Checksum

/**
 * Capability interface implemented by cloud-storage filesystem providers that
 * can return a content-addressable checksum for a path via a single metadata
 * call (e.g. S3 HEAD with {@code ChecksumMode.ENABLED}).
 *
 * Used by {@code FileHolder.funnel} when the active
 * {@link nextflow.util.CacheHelper.HashMode#CLOUD_HASH} mode is in effect, so
 * that task hashes contribute the cloud-native checksum (e.g. S3 CRC64NVMe)
 * instead of {@code path + size + mtime}. Filesystems that cannot provide a
 * checksum return {@link Optional#empty()} and the caller falls back to the
 * default hashing path.
 *
 * Mirrors the existing {@link FileSystemTransferAware} capability pattern:
 * checked at runtime via {@code instanceof} on the provider returned by
 * {@code path.getFileSystem().provider()}.
 */
interface ChecksumAwareFileSystemProvider {

    /**
     * Return the cloud-storage-native checksum for the given path, or
     * {@link Optional#empty()} when the path has no usable native checksum
     * (e.g. pre-2025 S3 object with null CRC64NVMe).
     *
     * Implementations MUST be cheap (a single HEAD-equivalent metadata call)
     * and side-effect-free.
     */
    Optional<Checksum> headChecksum(Path path)
}
