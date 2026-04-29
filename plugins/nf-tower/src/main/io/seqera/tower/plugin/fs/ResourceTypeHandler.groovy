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

package io.seqera.tower.plugin.fs

import java.nio.file.AccessMode
import java.nio.file.Path

/**
 * Strategy owning the semantics of one depth-3 path segment under {@code seqera://}.
 * Registered in {@link SeqeraFileSystem} at filesystem construction.
 *
 * Implementations own their resource's API client, caches, and interpretation of
 * trail segments beyond {@code resourceType}. The generic NIO layer does not look
 * inside the trail.
 */
interface ResourceTypeHandler {

    /** e.g. {@code "datasets"} or {@code "data-links"}. Must match the depth-3 path segment. */
    String getResourceType()

    /**
     * List entries at the given directory path. Caller has verified depth &ge; 3.
     * Returning an {@link Iterable} lets implementations stream large listings
     * without materializing them in memory.
     */
    Iterable<Path> list(SeqeraPath dir) throws IOException

    /** Return attributes for any path at depth &ge; 3 owned by this handler. */
    SeqeraFileAttributes readAttributes(SeqeraPath path) throws IOException

    /**
     * Open a read stream for a leaf path. Throw {@link IllegalArgumentException}
     * if the path is a directory or not otherwise addressable as a file.
     */
    InputStream newInputStream(SeqeraPath path) throws IOException

    /**
     * Verify the path exists and requested modes are satisfiable. READ is allowed;
     * WRITE/EXECUTE throw {@link java.nio.file.AccessDeniedException}.
     */
    void checkAccess(SeqeraPath path, AccessMode... modes) throws IOException
}
