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

import java.nio.file.DirectoryStream
import java.nio.file.Path

/**
 * Strategy interface for resource type handling within a Seqera workspace path.
 * Implementations provide listing and I/O for a specific resource type (e.g. "datasets").
 * Adding support for a new resource type requires only a new implementation of this interface.
 *
 * @author Seqera Labs
 */
interface ResourceTypeHandler {

    /**
     * @return the resource type string as it appears in the path (e.g. "datasets")
     */
    String getResourceType()

    /**
     * List entries under the given resource-type directory path.
     *
     * @param parent a depth-3 SeqeraPath representing the resource-type directory
     * @return a DirectoryStream of depth-4 paths (individual resource entries)
     */
    DirectoryStream<Path> listEntries(SeqeraPath parent) throws IOException

    /**
     * Open an InputStream for the given resource path.
     *
     * @param path a depth-4 SeqeraPath representing a specific resource
     * @return InputStream over the resource content
     */
    InputStream openInputStream(SeqeraPath path) throws IOException

    /**
     * Open an OutputStream for the given resource path.
     *
     * @param path a depth-4 SeqeraPath representing a target resource
     * @return OutputStream that uploads on close
     */
    OutputStream openOutputStream(SeqeraPath path) throws IOException
}
