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
 */

package io.seqera.tower.plugin.fs;

import nextflow.file.FileSystemPathFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Path;

/**
 * Factory for creating Seqera Platform Data-Link file system paths
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class SeqeraPathFactory extends FileSystemPathFactory {

    private static final Logger log = LoggerFactory.getLogger(SeqeraPathFactory.class);

    @Override
    protected Path parseUri(String uri) {
        if (!uri.startsWith(SeqeraPath.SEQERA_PROT)) {
            return null;
        }

        try {
            final SeqeraFileSystemProvider provider = new SeqeraFileSystemProvider();
            final java.net.URI parsedUri = SeqeraPath.asUri(uri);
            return provider.getPath(parsedUri);
        } catch (Exception e) {
            log.warn("Failed to parse Seqera URI: {}", uri, e);
            return null;
        }
    }

    @Override
    protected String toUriString(Path path) {
        if (path instanceof SeqeraPath) {
            return ((SeqeraPath) path).toUriString();
        }
        return null;
    }

    @Override
    protected String getBashLib(Path path) {
        return null;
    }

    @Override
    protected String getUploadCmd(String source, Path target) {
        // Could implement a custom upload command here
        return null;
    }
}
