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
/**
 * Implements a repository provider for GitHub service
 *
 * See https://gitlab.com/
 * See http://doc.gitlab.com/ee/api/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GitlabRepositoryProvider extends RepositoryProvider {

    GitlabRepositoryProvider(String project, ProviderConfig config=null) {
        this.project = project
        this.config = config ?: new ProviderConfig('gitlab')
    }

    final String getProjectName() {
        URLEncoder.encode(project,'utf-8')
    }

    protected void auth( HttpURLConnection connection ) {
        if( config.password ) {
            connection.setRequestProperty("PRIVATE-TOKEN", config.password)
        }
    }

    @Override
    String getName() { "GitLab" }

    @Override
    String getRepoUrl() {
        return "${config.endpoint}/api/v3/projects/${getProjectName()}"
    }

    String getDefaultBranch() {
        invokeAndParseResponse(getRepoUrl()) ?. default_branch
    }

    /** {@inheritDoc} */
    @Override
    String getContentUrl( String path ) {
        "${config.endpoint}/api/v3/projects/${getProjectName()}/repository/files?file_path=${path}&ref=${getDefaultBranch()}"
    }

    /** {@inheritDoc} */
    @Override
    String getCloneUrl() {
        Map response = invokeAndParseResponse( getRepoUrl() )

        def result = response.get('http_url_to_repo')
        if( !result )
            throw new IllegalStateException("Missing clone URL for: $project")

        return result
    }

    /** {@inheritDoc} */
    @Override
    String getHomePage() {
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
