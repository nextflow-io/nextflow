/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
