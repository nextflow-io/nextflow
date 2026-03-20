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

import groovy.transform.CompileStatic
import io.seqera.tower.model.DatasetDto

/**
 * {@link ResourceTypeHandler} implementation for Seqera Platform datasets.
 * Lists datasets in a workspace and delegates I/O to {@link SeqeraFileSystemProvider}.
 *
 * @author Seqera Labs
 */
@CompileStatic
class DatasetsResourceHandler implements ResourceTypeHandler {

    @Override
    String getResourceType() { 'datasets' }

    @Override
    DirectoryStream<Path> listEntries(SeqeraPath parent) throws IOException {
        final fs = parent.getFileSystem() as SeqeraFileSystem
        final workspaceId = fs.resolveWorkspaceId(parent.org, parent.workspace)
        final datasets = fs.resolveDatasets(workspaceId)
        final List<Path> entries = datasets.collect { DatasetDto ds -> parent.resolve(ds.name) as Path }
        return new DirectoryStream<Path>() {
            @Override Iterator<Path> iterator() { entries.iterator() }
            @Override void close() {}
        }
    }

    @Override
    InputStream openInputStream(SeqeraPath path) throws IOException {
        return ((SeqeraFileSystemProvider) path.getFileSystem().provider()).newInputStream(path)
    }

    @Override
    OutputStream openOutputStream(SeqeraPath path) throws IOException {
        return ((SeqeraFileSystemProvider) path.getFileSystem().provider()).newOutputStream(path)
    }
}
