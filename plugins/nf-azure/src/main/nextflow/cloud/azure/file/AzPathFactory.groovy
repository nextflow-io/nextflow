/*
 * Copyright 2021, Microsoft Corp
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
import nextflow.cloud.azure.batch.AzFileCopyStrategy
import nextflow.cloud.azure.config.AzConfig
import nextflow.cloud.azure.nio.AzPath
import nextflow.file.FileHelper
import nextflow.file.FileSystemPathFactory
/**
 * Create Azure path objects for az:// prefixed URIs
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AzPathFactory extends FileSystemPathFactory {

    AzPathFactory() {
      log.debug "Creating Azure path factory"
    }

    @Override
    protected Path parseUri(String uri) {
        if( !uri.startsWith('az://') )
            return null

        if( uri.startsWith('az:///') )
            throw new IllegalArgumentException("Invalid Azure path URI - make sure the schema prefix does not container more than two slash characters - offending value: $uri")

        final storageConfigEnv = AzConfig.getConfig().storage().getEnv()
        final activeDirectoryConfigEnv = AzConfig.getConfig().activeDirectory().getEnv()
        final managedIdentityConfigEnv = AzConfig.getConfig().managedIdentity().getEnv()

        final configEnv = storageConfigEnv + activeDirectoryConfigEnv + managedIdentityConfigEnv

        // find the related file system
        final fs = getFileSystem(uri0(uri), configEnv)

        // resulting az path
        return fs.getPath(uri.substring(4))
    }

    private URI uri0(String uri) {
        // note: this is needed to allow URI to handle curly brackets characters
        // see https://github.com/nextflow-io/nextflow/issues/1969
        new URI(null, null, uri, null, null)
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
            log.trace "+ Setting context class loader to=$loader - previous=$bak"
            Thread.currentThread().setContextClassLoader(loader)
            return FileHelper.getOrCreateFileSystemFor(uri, env)
        }
        finally {
            Thread.currentThread().setContextClassLoader(bak)
        }
    }

    @Override
    protected String toUriString(Path path) {
        return path instanceof AzPath ? ((AzPath)path).toUriString() : null
    }

    @Override
    protected String getBashLib(Path path) {
        return path instanceof AzPath ? AzBashLib.script() : null
    }

    @Override
    protected String getUploadCmd(String source, Path target) {
        return target instanceof AzPath ?  AzFileCopyStrategy.uploadCmd(source, target) : null
    }

}
