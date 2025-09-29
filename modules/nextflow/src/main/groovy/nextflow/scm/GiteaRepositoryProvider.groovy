/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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


import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 * Implements a repository provider for Gitea service
 *
 * @author Akira Sekiguchi <pachiras.yokohama@gmail.com>
 */
@Slf4j
@CompileStatic
final class GiteaRepositoryProvider extends RepositoryProvider {

    GiteaRepositoryProvider(String project, ProviderConfig config=null) {
        this.project = project
        this.config = config ?: new ProviderConfig('gitea')
    }

    /** {@inheritDoc} */
    @Override
    String getName() { "Gitea" }

    /** {@inheritDoc} */
    @Override
    String getEndpointUrl() {
        "${config.endpoint}/repos/${project}"
    }

    @Override
    protected String[] getAuth() {
        if( config.token ) {
            // set the token in the request header
            // https://docs.gitea.io/en-us/api-usage/#authentication
            return new String[] { "Authorization", "token $config.token" as String }
        }
        else
            return null
    }

    @Override
    boolean hasCredentials() {
        return getToken()
            ? true
            : super.hasCredentials()
    }

    @Override
    @CompileDynamic
    List<BranchInfo> getBranches() {
        // https://try.gitea.io/api/swagger#/repository/repoListBranches
        final url = "${config.endpoint}/repos/${project}/branches"
        this.<BranchInfo>invokeAndResponseWithPaging(url, { Map branch -> new BranchInfo(branch.name as String, branch.commit?.id as String) })
    }

    @Override
    @CompileDynamic
    List<TagInfo> getTags() {
        // https://try.gitea.io/api/swagger#/repository/repoListTags
        final url = "${config.endpoint}/repos/${project}/tags"
        this.<TagInfo>invokeAndResponseWithPaging(url, { Map tag -> new TagInfo(tag.name as String, tag.commit?.sha as String) })
    }

    /** {@inheritDoc} */
    @Override
    String getContentUrl( String path ) {
        // see
        // https://try.gitea.io/api/swagger#/repository/repoGetRawFile
        // note: `ref` is undocumented
        def result = "${config.endpoint}/repos/$project/raw/$path"
        if( revision )
            result += "?ref=$revision"
        return result
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
        return invokeBytes(url)
    }

    /** {@inheritDoc} */
    @Override
    List<RepositoryEntry> listDirectory(String path, int depth) {
        final branch = revision ?: "master"
        // Normalize path using base class helper
        final dirPath = normalizePath(path)
        
        // Build the contents API URL - Gitea follows GitHub-like API pattern
        String url = "${config.endpoint}/repos/$project/contents"
        if (dirPath) {
            url += "/$dirPath"
        }
        url += "?ref=$branch"
        
        try {
            // Make the API call
            def response = invoke(url)
            List<Map> contents = new groovy.json.JsonSlurper().parseText(response) as List<Map>
            
            if (!contents) {
                return []
            }
            
            List<RepositoryEntry> entries = []
            
            for (Map entry : contents) {
                // For depth control, we need to filter appropriately
                if (shouldIncludeEntry(entry, path, depth)) {
                    entries.add(createRepositoryEntry(entry))
                }
            }
            
            // If depth > 1, we need to recursively get subdirectory contents
            if (depth > 1) {
                entries.addAll(getRecursiveEntries(path, depth, branch, 1))
            }
            
            return entries.sort { it.name }
            
        } catch (Exception e) {
            throw new UnsupportedOperationException("Directory listing failed for Gitea path: $path", e)
        }
    }

    private List<RepositoryEntry> getRecursiveEntries(String basePath, int maxDepth, String branch, int currentDepth) {
        if (currentDepth > maxDepth) {
            return []
        }
        
        List<RepositoryEntry> allEntries = []
        
        // Get current level entries first
        final normalizedBasePath = normalizePath(basePath)
        String url = "${config.endpoint}/repos/$project/contents"
        if (normalizedBasePath) {
            url += "/$normalizedBasePath"
        }
        url += "?ref=$branch"
        
        try {
            def response = invoke(url)
            List<Map> contents = new groovy.json.JsonSlurper().parseText(response) as List<Map>
            
            for (Map entry : contents) {
                if (entry.get('type') == 'dir' && currentDepth < maxDepth) {
                    String entryName = entry.get('name') as String
                    String subPath = normalizedBasePath ? "$normalizedBasePath/$entryName" : entryName
                    allEntries.addAll(getRecursiveEntries(subPath, maxDepth, branch, currentDepth + 1))
                }
            }
        } catch (Exception e) {
            log.debug("Failed to process directory during recursive listing: ${e.message}")
            // Continue processing other directories if one fails
        }
        
        return allEntries
    }

    private boolean shouldIncludeEntry(Map entry, String basePath, int depth) {
        // For Gitea contents API, we get direct children, so depth 0 includes all returned items
        return true
    }

    private RepositoryEntry createRepositoryEntry(Map entry) {
        String name = entry.get('name') as String
        String path = entry.get('path') as String
        String type = entry.get('type') as String
        
        EntryType entryType = (type == 'dir') ? EntryType.DIRECTORY : EntryType.FILE
        String sha = entry.get('sha') as String
        Long size = entry.get('size') as Long
        
        // Ensure absolute path using base class helper
        String fullPath = ensureAbsolutePath(path)
        
        return new RepositoryEntry(
            name: name,
            path: fullPath,
            type: entryType,
            sha: sha,
            size: size
        )
    }

}
