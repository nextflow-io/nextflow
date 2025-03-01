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
 *
 */

package nextflow.data.cid.fs

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.data.config.DataConfig
import nextflow.file.FileHelper
import nextflow.file.FileSystemPathFactory

import static nextflow.data.cid.fs.CidPath.CID_PROT

/**
 * Implements a {@link FileSystemPathFactory} for CID file system
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
class CidPathFactory extends FileSystemPathFactory {

    @Override
    protected Path parseUri(String uri) {
        return uri.startsWith(CID_PROT) ? create(uri) : null
    }

    @Override
    protected String toUriString(Path path) {
        return path instanceof CidPath ? ((CidPath)path).toUriString() : null
    }

    @Override
    protected String getBashLib(Path target) {
        return null
    }

    @Override
    protected String getUploadCmd(String source, Path target) {
        return null
    }

    static CidPath create(String path) {
        final uri = CidPath.asUri(path)
        return (CidPath) FileHelper.getOrCreateFileSystemFor(uri, DataConfig.asMap()).provider().getPath(uri)
    }
}
