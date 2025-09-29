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
 */

package nextflow.scm

import groovy.json.JsonSlurper
import groovy.util.logging.Slf4j
/**
 * Implements a repository provider for GitHub service
 *
 * See https://gitlab.com/
 * https://docs.gitlab.com/ee/api/repositories.html
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class GitlabRepositoryProvider extends RepositoryProvider {

    GitlabRepositoryProvider(String project, ProviderConfig config=null) {
        this.project = project
        this.config = config ?: new ProviderConfig('gitlab')
    }

    final String getProjectName() {
        URLEncoder.encode(project,'utf-8')
    }

    @Override
    protected String[] getAuth() {
        if( config.token ) {
            // set the token in the request header
            return new String[] { "PRIVATE-TOKEN", config.token }
        }
        if( config.password ) {
            return new String[] { "PRIVATE-TOKEN", config.password }
        }
        return null
    }

    @Override
    boolean hasCredentials() {
        return getToken()
            ? true
            : super.hasCredentials()
    }

    @Override
    String getName() { "GitLab" }

    @Override
    String getEndpointUrl() {
        return "${config.endpoint}/api/v4/projects/${getProjectName()}"
    }

    String getDefaultBranch() {
        def result = invokeAndParseResponse(getEndpointUrl()) ?. default_branch
        if( !result ) {
            log.debug "Unable to fetch repo default branch. Using `master` branch -- See https://gitlab.com/gitlab-com/support-forum/issues/1655#note_26132691"
            return 'master'
        }
        return result
    }

    @Override
    List<BranchInfo> getBranches() {
        // https://docs.gitlab.com/ee/api/branches.html
        final url = "${config.endpoint}/api/v4/projects/${getProjectName()}/repository/branches"
        this.<BranchInfo>invokeAndResponseWithPaging(url, { Map branch -> new BranchInfo(branch.name as String, branch.commit?.id as String) })
    }

    @Override
    List<TagInfo> getTags() {
        // https://docs.gitlab.com/ee/api/tags.html
        final url = "${config.endpoint}/api/v4/projects/${getProjectName()}/repository/tags"
        this.<TagInfo>invokeAndResponseWithPaging(url, { Map tag -> new TagInfo(tag.name as String, tag.commit?.id as String) })
    }

    /** {@inheritDoc} */
    @Override
    String getContentUrl( String path ) {
        // see
        //  https://docs.gitlab.com/ee/api/repository_files.html#get-raw-file-from-repository
        //
        final ref = revision ?: getDefaultBranch()
        final encodedPath = URLEncoder.encode(path.stripStart('/'),'utf-8')
        return "${config.endpoint}/api/v4/projects/${getProjectName()}/repository/files/${encodedPath}?ref=${ref}"
    }

    /** {@inheritDoc} */
    @Override
    String getCloneUrl() {
        Map response = invokeAndParseResponse( getEndpointUrl() )

        def result = response.get('http_url_to_repo')
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

    /** {@inheritDoc} */
    @Override
    List<RepositoryEntry> listDirectory(String path, int depth) {
        final ref = revision ?: getDefaultBranch()
        final encodedPath = path ? URLEncoder.encode(path.stripStart('/'), 'utf-8') : ""
        
        // Build the Tree API URL
        String url = "${config.endpoint}/api/v4/projects/${getProjectName()}/repository/tree"
        List<String> params = []
        if (ref) params.add("ref=${ref}")
        if (encodedPath) params.add("path=${encodedPath}")
        
        // For GitLab, we use recursive=true for any depth > 1
        if (depth > 1) {
            params.add("recursive=true")
        }
        
        if (params) {
            url += "?" + params.join("&")
        }
        
        // Make the API call and parse response
        String response = invoke(url)
        List<Map> treeEntries = response ? new JsonSlurper().parseText(response) as List<Map> : []
        
        if (!treeEntries) {
            return []
        }
        
        List<RepositoryEntry> entries = []
        
        for (Map entry : treeEntries) {
            String entryPath = entry.get('path') as String
            
            // Filter entries based on depth
            if (shouldIncludeEntry(entryPath, path, depth)) {
                entries.add(createRepositoryEntry(entry, path))
            }
        }
        
        return entries.sort { it.name }
    }

    private boolean shouldIncludeEntry(String entryPath, String basePath, int depth) {
        String relativePath = entryPath
        if (basePath && !basePath.isEmpty()) {
            // If we have a base path, compute the relative path
            String normalizedBase = basePath.stripStart('/').stripEnd('/')
            String normalizedEntry = entryPath.stripStart('/').stripEnd('/')
            
            if (normalizedEntry.startsWith(normalizedBase + "/")) {
                relativePath = normalizedEntry.substring(normalizedBase.length() + 1)
            } else if (normalizedEntry == normalizedBase) {
                return false // Skip the base directory itself
            }
        }
        
        // Count directory levels in the relative path
        int entryDepth = relativePath.split("/").length - 1
        
        // Include if within depth limit: depth=1 includes immediate children only,
        // depth=2 includes children+grandchildren, depth=3 includes children+grandchildren+great-grandchildren, etc.
        return entryDepth < depth
    }

    private RepositoryEntry createRepositoryEntry(Map entry, String basePath) {
        String entryPath = entry.get('path') as String
        String name = entry.get('name') as String
        
        EntryType type = entry.get('type') == 'tree' ? EntryType.DIRECTORY : EntryType.FILE
        String sha = entry.get('id') as String
        Long size = null // GitLab tree API doesn't provide file size
        
        return new RepositoryEntry(
            name: name,
            path: entryPath,
            type: type,
            sha: sha,
            size: size
        )
    }
}
