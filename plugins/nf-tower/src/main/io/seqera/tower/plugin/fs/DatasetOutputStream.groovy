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

import java.nio.file.NoSuchFileException

import groovy.transform.CompileStatic

/**
 * Buffered {@link OutputStream} that uploads to the Seqera Platform on {@code close()}.
 * Enforces the 10 MB platform size limit before making the API call.
 *
 * @author Seqera Labs
 */
@CompileStatic
class DatasetOutputStream extends OutputStream {

    static final long MAX_UPLOAD_BYTES = 10L * 1024 * 1024  // 10 MB

    private final SeqeraPath path
    private final ByteArrayOutputStream buffer = new ByteArrayOutputStream()
    private boolean closed = false

    DatasetOutputStream(SeqeraPath path) {
        this.path = path
    }

    @Override
    void write(int b) throws IOException {
        buffer.write(b)
    }

    @Override
    void write(byte[] b, int off, int len) throws IOException {
        buffer.write(b, off, len)
    }

    @Override
    void close() throws IOException {
        if (closed) return
        closed = true
        final bytes = buffer.toByteArray()
        if (bytes.length > MAX_UPLOAD_BYTES)
            throw new IOException("Dataset '${path.datasetName}' exceeds the 10 MB Seqera Platform limit (${bytes.length} bytes)")

        final fs = path.getFileSystem() as SeqeraFileSystem
        final workspaceId = fs.resolveWorkspaceId(path.org, path.workspace)

        // Resolve or create the dataset
        String datasetId
        try {
            final existing = fs.resolveDataset(workspaceId, path.datasetName)
            datasetId = existing.id
        } catch (NoSuchFileException ignored) {
            final created = fs.client.createDataset(workspaceId, path.datasetName)
            datasetId = created.id
        }

        final fileName = path.datasetName
        fs.client.uploadDataset(datasetId, bytes, fileName, false)
        fs.invalidateDatasetCache(workspaceId)
    }
}
