/*
 * Copyright 2013-2024, Seqera Labs
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

package io.seqera.tower.plugin.dataset

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.FileHelper
import nextflow.file.FileSystemPathFactory

/**
 * Factory that intercepts {@code dataset://} URI strings and returns
 * a {@link DatasetPath} via the NIO FileSystem provider.
 * <p>
 * Registered as a Nextflow {@code ExtensionPoint} so that
 * {@code FileHelper.asPath("dataset://my-samplesheet")} works
 * transparently — no pipeline code changes needed.
 *
 * @author Edmund Miller
 */
@Slf4j
@CompileStatic
class DatasetPathFactory extends FileSystemPathFactory {

    @Override
    protected Path parseUri(String str) {
        if (!str.startsWith('dataset://'))
            return null

        log.debug "Parsing dataset URI: {}", str

        // Normalise to triple-slash form for URI parsing:
        // dataset://name → dataset:///name
        final normalized = str.startsWith('dataset:///') ? str : 'dataset:///' + str.substring('dataset://'.length())

        final uri = new URI(null, null, normalized, null, null)
        return FileHelper.getOrCreateFileSystemFor(uri).provider().getPath(uri)
    }

    @Override
    protected String toUriString(Path path) {
        if (path instanceof DatasetPath) {
            return path.toString()
        }
        return null
    }

    @Override
    protected String getBashLib(Path target) {
        // dataset:// paths are resolved to cloud paths before execution,
        // no special bash lib needed
        return null
    }

    @Override
    protected String getUploadCmd(String source, Path target) {
        // read-only — no upload support
        return null
    }
}
