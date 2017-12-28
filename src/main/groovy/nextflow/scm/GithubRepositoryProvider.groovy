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

import groovy.transform.CompileStatic

/**
 * Implements a repository provider for GitHub service
 *
 * @author Maria Chatzou
 * @author Paolo Di Tommaso
 */
@CompileStatic
final class GithubRepositoryProvider extends RepositoryProvider {

    GithubRepositoryProvider(String project, ProviderConfig config=null) {
        this.project = project
        this.config = config ?: new ProviderConfig('github')
    }

    /** {@inheritDoc} */
    @Override
    String getName() { "GitHub" }

    /** {@inheritDoc} */
    @Override
    String getEndpointUrl() {
        "${config.endpoint}/repos/${project}"
    }

    /** {@inheritDoc} */
    @Override
    String getContentUrl( String path ) {
        "${config.endpoint}/repos/$project/contents/$path"
    }

    /** {@inheritDoc} */
    @Override
    String getCloneUrl() {
        Map response = invokeAndParseResponse( getEndpointUrl() )

        def result = response.get('clone_url')
        if( !result )
            throw new IllegalStateException("Missing clone URL for: $project")

        return result
    }

    /** {@inheritDoc} */
    @Override
    String getRepositoryUrl() {
        "${config.server}/$project"
    }

    /** {@inheritDoc} */
    @Override
    byte[] readBytes(String path) {

        def url = getContentUrl(path)
        Map response  = invokeAndParseResponse(url)
        response.get('content')?.toString()?.decodeBase64()

    }

}
