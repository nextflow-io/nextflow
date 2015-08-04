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

import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
/**
 * Implements a repository provider for the BitBucket service
 *
 * Author Maria Chatzou
 * Author Paolo Di Tommaso
 */
@Slf4j
final class BitbucketRepositoryProvider extends RepositoryProvider {

    /** {@inheritDoc} */
    @Override
    String getName() { "BitBucket" }

    @Override
    String getRepoUrl() {
        return "https://bitbucket.org/api/2.0/repositories/${pipeline}"
    }

    @Override
    String getContentUrl( String path ) {
        return "https://bitbucket.org/api/1.0/repositories/$pipeline/src/${getMainBranch()}/$path"
    }

    private String getMainBranchUrl() {
        "https://bitbucket.org/api/1.0/repositories/$pipeline/main-branch/"
    }

    String getMainBranch() {
        invokeAndParseResponse(getMainBranchUrl()) ?. name
    }

    @Override
    String getCloneUrl() {
        Map response = invokeAndParseResponse( getRepoUrl() )

        if( response?.scm != "git" ){
            throw new AbortOperationException("Bitbucket repository at ${getHomePage()} is not supporting Git")
        }

        def result = response?.links?.clone?.find{ it.name == "https" } as Map
        if( !result )
            throw new IllegalStateException("Missing clone URL for: $pipeline")

        return result.href
    }

    @Override
    String getHomePage() {
        return "https://bitbucket.org/$pipeline"
    }

    @Override
    byte[] readBytes(String path) {

        def url = getContentUrl(path)
        Map response  = invokeAndParseResponse(url)
        response.get('data')?.toString()?.getBytes()

    }
}