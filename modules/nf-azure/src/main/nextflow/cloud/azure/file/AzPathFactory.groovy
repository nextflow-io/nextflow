/*
 * Copyright 2020, Microsoft Corp
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

package nextflow.cloud.azure.file

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.cloud.azure.config.AzConfig
import nextflow.file.FileHelper
import nextflow.file.FileSystemPathFactory
/**
 * Create Azure path objects for azb:// prefixed URIs
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AzPathFactory extends FileSystemPathFactory {

    @Override
    protected Path parseUri(String uri) {
        if( !uri.startsWith('azb://') )
            return null

        final cfg = AzConfig.getConfig().storage()

        // parse uri
        final tokens = AzContainerTokens.parse(uri)
        final account = tokens.account ?: cfg.accountName
        if( !tokens.bucket )
            throw new IllegalArgumentException("Invalid Azure storage container URI: $uri")

        // find the related file system 
        final String accountUri = "azb://?account=$account"
        final fs = FileHelper.getOrCreateFileSystemFor(new URI(accountUri), cfg.getEnv())

        // compose the target path
        return tokens.path ? fs.getPath(tokens.bucket, tokens.path) : fs.getPath(tokens.bucket)
    }


    @Override
    protected String toUriString(Path path) {
        return null
    }


}
