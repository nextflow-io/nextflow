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

import java.nio.file.FileSystem
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.azure.AzurePlugin
import nextflow.cloud.azure.config.AzConfig
import nextflow.file.FileHelper
import nextflow.file.FileSystemPathFactory
/**
 * Create Azure path objects for azb:// prefixed URIs
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AzPathFactory extends FileSystemPathFactory {

    AzPathFactory() {
      log.debug "Creating Azure path fatory"
    }

    @Override
    protected Path parseUri(String uri) {
        if( !uri.startsWith('azb://') )
            return null

        final cfg = AzConfig.getConfig().storage()

        // parse uri
        final tokens = AzStorageContainerParser.parse(uri)
        final account = tokens.account ?: cfg.accountName
        if( !tokens.container )
            throw new IllegalArgumentException("Invalid Azure storage container URI: $uri")

        // find the related file system 
        final String accountUri = "azb://?account=$account"
        final fs = getFileSystem(new URI(accountUri), cfg.getEnv())

        // compose the target path
        return tokens.path ? fs.getPath(tokens.container+':', tokens.path) : fs.getPath(tokens.container)
    }

    protected FileSystem getFileSystem(URI uri, Map env) {
        final bak = Thread.currentThread().getContextClassLoader()
        // NOTE: setting the context classloader to allow loading azure deps via java ServiceLoader
        // see
        //  com.azure.core.http.HttpClientProvider
        //  com.azure.core.http.netty.implementation.ReactorNettyClientProvider
        //
        try {
            final loader = AzurePlugin.class.getClassLoader()
            log.debug "+ Setting context class loader to=$loader - previous=$bak"
            Thread.currentThread().setContextClassLoader(loader)
            return FileHelper.getOrCreateFileSystemFor(uri, env)
        }
        finally {
            Thread.currentThread().setContextClassLoader(bak)
        }
    }

    @Override
    protected String toUriString(Path path) {
        return 'azb://' + path.toAbsolutePath().toString()
    }

}
