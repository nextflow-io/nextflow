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

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.SysEnv

import java.nio.charset.StandardCharsets

/**
 * Implements a repository provider for GitHub service
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class GithubRepositoryProvider extends RepositoryProvider {

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
        return "${config.endpoint}/repos/${project}"
    }

    @Override
    boolean hasCredentials() {
        super.hasCredentials() ?: SysEnv.containsKey('GITHUB_TOKEN')
    }

    @Override
    String getUser() {
        super.getUser() ?: SysEnv.get('GITHUB_TOKEN')
    }

    @Override
    String getPassword() {
        super.getPassword() ?: (SysEnv.containsKey('GITHUB_TOKEN') ? 'x-oauth-basic' : null)
    }

    /** {@inheritDoc} */
    @Override
    String getContentUrl( String path ) {
        // see
        // https://docs.github.com/en/rest/reference/repos#get-repository-content
        //
        def result = "${config.endpoint}/repos/$project/contents/$path"
        if( revision )
            result += "?ref=${URLEncoder.encode(revision, StandardCharsets.UTF_8)}"
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

    @Override
    @CompileDynamic
    @Memoized
    List<BranchInfo> getBranches() {
        // https://developer.github.com/v3/repos/branches/#list-branches
        final url = "${config.endpoint}/repos/$project/branches"
        this.<BranchInfo>invokeAndResponseWithPaging(url, { Map branch -> new BranchInfo(branch.name as String, branch.commit?.sha as String) })
    }

    @Override
    @CompileDynamic
    @Memoized
    List<TagInfo> getTags() {
        // https://developer.github.com/v3/repos/#list-tags
        final url = "${config.endpoint}/repos/$project/tags"
        this.<TagInfo>invokeAndResponseWithPaging(url, { Map tag -> new TagInfo(tag.name as String, tag.commit?.sha as String)})
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
    @Memoized
    List<RepositoryEntry> listDirectory(String path, int depth) {
        // Get the tree SHA for the specific directory
        String treeSha = getTreeSha(path)
        
        // Build the Trees API URL
        String url = getTreeUrl(treeSha, depth > 1)
        
        // Make the API call and parse response
        Map response = invokeAndParseResponse(url)
        List<Map> treeEntries = response.get('tree') as List<Map>
        
        if (!treeEntries) {
            return []
        }
        
        List<RepositoryEntry> entries = []
        
        for (Map entry : treeEntries) {
            String entryPath = entry.get('path') as String
            
            // Include if within depth limit: depth=0 includes immediate children only,
            // depth=1 includes children+grandchildren, depth=2 includes children+grandchildren+great-grandchildren, etc.
            int entryDepth = entryPath.split("/").length - 1
            if (depth == -1 || entryDepth <= depth) {
                entries.add(createRepositoryEntry(entry, path))
            }
        }
        
        return entries.sort { it.name }
    }

    private String getTreeUrl(String treeSha, boolean recursive) {
        String url = "${config.endpoint}/repos/$project/git/trees/$treeSha"
        if (recursive) {
            url += "?recursive=1"
        }
        return url
    }

    @Memoized
    private String getTreeSha(String path) {
        // Normalize path using base class helper
        def normalizedPath = normalizePath(path)
        
        if (normalizedPath && !normalizedPath.isEmpty()) {
            // For subdirectory, we need to find the tree SHA by traversing from root
            return getTreeShaForPath(normalizedPath)
        }
        
        // For root directory, get the commit SHA and then the tree SHA
        String commitSha = getCommitSha()
        Map commit = invokeAndParseResponse("${config.endpoint}/repos/$project/git/commits/$commitSha")
        Map tree = commit.get('tree') as Map
        return tree.get('sha') as String
    }

    private String getTreeShaForPath(String path) {
        // Start from root tree
        String currentTreeSha = getTreeSha("")
        String[] pathParts = path.split("/")
        
        for (String part : pathParts) {
            String url = getTreeUrl(currentTreeSha, false)
            Map response = invokeAndParseResponse(url)
            List<Map> treeEntries = response.get('tree') as List<Map>
            
            Map foundEntry = treeEntries.find { 
                it.get('path') == part && it.get('type') == 'tree'
            }
            
            if (!foundEntry) {
                throw new IllegalArgumentException("Directory not found: $path")
            }
            
            currentTreeSha = foundEntry.get('sha') as String
        }
        
        return currentTreeSha
    }

    @Memoized
    private String getCommitSha() {
        if (revision) {
            // Try to resolve the revision to a commit SHA
            try {

                Map ref = invokeAndParseResponse("${config.endpoint}/repos/$project/git/refs/heads/${URLEncoder.encode(revision, StandardCharsets.UTF_8)}")
                Map object = ref.get('object') as Map
                return object.get('sha') as String
            } catch (Exception e) {
                // If it's not a branch, try as a tag or direct SHA
                return revision
            }
        }
        
        // Default to main/master branch
        try {
            Map ref = invokeAndParseResponse("${config.endpoint}/repos/$project/git/refs/heads/main")
            Map object = ref.get('object') as Map
            return object.get('sha') as String
        } catch (Exception e) {
            Map ref = invokeAndParseResponse("${config.endpoint}/repos/$project/git/refs/heads/master")
            Map object = ref.get('object') as Map
            return object.get('sha') as String
        }
    }


    private RepositoryEntry createRepositoryEntry(Map entry, String basePath) {
        String entryPath = entry.get('path') as String
        
        // Create absolute path using base class helper
        def normalizedBasePath = normalizePath(basePath)
        String fullPath = normalizedBasePath && !normalizedBasePath.isEmpty() ? "/${normalizedBasePath}/${entryPath}" : ensureAbsolutePath(entryPath)
        
        // For name, use just the entry path (which is relative to the directory we're listing)
        String name = entryPath
        if (entryPath.contains("/")) {
            name = entryPath.substring(entryPath.lastIndexOf("/") + 1)
        }
        
        EntryType type = entry.get('type') == 'tree' ? EntryType.DIRECTORY : EntryType.FILE
        String sha = entry.get('sha') as String
        Long size = entry.get('size') as Long
        
        return new RepositoryEntry(
            name: name,
            path: fullPath,
            type: type,
            sha: sha,
            size: size
        )
    }

}
