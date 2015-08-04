/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

import groovy.transform.CompileStatic

/**
 * Implements a repository provider for GitHub service
 *
 * Author Maria Chatzou
 * Author Paolo Di Tommaso
 */
@CompileStatic
final class GithubRepositoryProvider extends RepositoryProvider{

    /** {@inheritDoc} */
    @Override
    String getName() { "GitHub" }

    /** {@inheritDoc} */
    @Override
    String getRepoUrl() {
        "https://api.github.com/repos/${pipeline}"
    }

    /** {@inheritDoc} */
    @Override
    String getContentUrl( String path ) {
        "https://api.github.com/repos/$pipeline/contents/$path"
    }

    /** {@inheritDoc} */
    @Override
    String getCloneUrl() {
        Map response = invokeAndParseResponse( getRepoUrl() )

        def result = response.get('clone_url')
        if( !result )
            throw new IllegalStateException("Missing clone URL for: $pipeline")

        return result
    }

    /** {@inheritDoc} */
    @Override
    String getHomePage() {
        "https://github.com/$pipeline"
    }

    /** {@inheritDoc} */
    @Override
    byte[] readBytes(String path) {

        def url = getContentUrl(path)
        Map response  = invokeAndParseResponse(url)
        response.get('content')?.toString()?.decodeBase64()

    }
}
