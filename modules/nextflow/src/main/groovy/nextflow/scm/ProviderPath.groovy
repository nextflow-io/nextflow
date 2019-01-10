/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.scm

import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.Memoized
import groovy.transform.PackageScope

/**
 * Implements a {@link Path} interface for the a given
 * source code {@link RepositoryProvider}
 *
 * NOTE: this is not suppose to full support NIO2 {@link java.nio.file.Files} API
 * It designed only to allow the reading of a config files included in a remote
 * project hosted in a source code repository
 *
 * @see nextflow.config.ConfigParser
 * @see nextflow.config.ConfigBase
 *
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProviderPath implements Path {

    @PackageScope
    RepositoryProvider provider

    @PackageScope
    @Delegate
    Path delegate

    ProviderPath(RepositoryProvider provider, String path) {
        this(provider, Paths.get(path))
    }

    private ProviderPath(RepositoryProvider provider, Path path) {
        assert provider
        assert !path.isAbsolute()
        this.provider = provider
        this.delegate = path
    }

    Path resolveSibling(Path other) {
        new ProviderPath(provider, delegate.resolveSibling(other.toString()))
    }

    Path resolveSibling(String other) {
        new ProviderPath(provider, delegate.resolveSibling(other))
    }

    Path resolve(Path other) {
        new ProviderPath(provider, delegate.resolve(other.toString()))
    }

    Path resolve(String other) {
        new ProviderPath(provider, delegate.resolve(other))
    }

    boolean exists(LinkOption... options) { true }

    @Memoized
    String getText() {
        provider.readText(delegate.toString())
    }

    String toUriString() {
        provider.getContentUrl(delegate.toString())
    }

}
