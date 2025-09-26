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

import org.eclipse.jgit.api.Git
import org.eclipse.jgit.lib.Constants
import org.eclipse.jgit.lib.Ref
import org.eclipse.jgit.revwalk.RevWalk
import org.eclipse.jgit.treewalk.TreeWalk
/**
 * Local storage asset provider
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LocalRepositoryProvider extends RepositoryProvider {

    LocalRepositoryProvider( String project, ProviderConfig config) {
        assert project
        assert config
        assert config.path
        assert !project.contains('/')

        this.project = project
        this.config = config
        this.path = new File(config.path)
    }

    private File path

    @Override
    String getName() { 'local' }

    @Override
    String getEndpointUrl() {
        return "file:${new File(path, project)}"
    }

    @Override
    String getContentUrl(String path) {
        def base = new File(this.path, project)
        def result = new File(base, path)
        return "file:$result"
    }

    @Override
    String getCloneUrl() {
        final root = new File(path, project)
        return new File(root,'.git').isDirectory() ? "file:${root}/.git" : "file:${root}"
    }

    @Override
    String getRepositoryUrl() {
        new File(path, project).toString()
    }


    @Override
    byte[] readBytes(String path) {

        final git = Git.open(new File(this.path, project))
        try {
            final repo = git.getRepository()

            def lastCommitId = repo.resolve(Constants.HEAD)
            def revWalk = new RevWalk(repo)
            def commit = revWalk.parseCommit(lastCommitId)
            def tree = commit.getTree()

            def treeWalk = TreeWalk.forPath(repo, path, tree)
            if( !treeWalk )
                return null
            def id = treeWalk.getObjectId(0)
            def loader = repo.open(id)

            def source = loader.openStream()
            def result = new ByteArrayOutputStream()
            int ch
            while( (ch=source.read()) != -1 ) {
                result.write(ch)
            }
            result.close()
            source.close()
            treeWalk.close()

            return result.toByteArray()
        }
        finally {
            git.close()
        }
    }

    @Override
    List<RepositoryEntry> listDirectory(String path, int depth) {
        final git = Git.open(new File(this.path, project))
        try {
            final repo = git.getRepository()
            def lastCommitId = repo.resolve(Constants.HEAD)
            def revWalk = new RevWalk(repo)
            def commit = revWalk.parseCommit(lastCommitId)
            def tree = commit.getTree()
            
            def treeWalk = new TreeWalk(repo)
            treeWalk.addTree(tree)
            treeWalk.setRecursive(depth != 0)
            
            // Set path filter if specified
            if (path && !path.isEmpty()) {
                treeWalk.setFilter(org.eclipse.jgit.treewalk.filter.PathFilter.create(path))
            }
            
            List<RepositoryEntry> entries = []
            
            while (treeWalk.next()) {
                String entryPath = treeWalk.getPathString()
                
                // Skip the base path itself
                if (path && entryPath == path) {
                    continue
                }
                
                // Filter by depth
                if (shouldIncludeEntry(entryPath, path, depth)) {
                    entries.add(createRepositoryEntry(treeWalk, entryPath))
                }
            }
            
            treeWalk.close()
            revWalk.close()
            
            return entries.sort { it.name }
            
        } finally {
            git.close()
        }
    }

    private boolean shouldIncludeEntry(String entryPath, String basePath, int depth) {
        if (depth == -1) {
            return true // Include all entries for fully recursive
        }
        
        String relativePath = entryPath
        if (basePath && !basePath.isEmpty()) {
            String normalizedBase = basePath.stripStart('/').stripEnd('/')
            String normalizedEntry = entryPath.stripStart('/').stripEnd('/')
            
            if (normalizedEntry.startsWith(normalizedBase + "/")) {
                relativePath = normalizedEntry.substring(normalizedBase.length() + 1)
            } else if (normalizedEntry == normalizedBase) {
                return false // Skip the base directory itself
            } else {
                return false // Entry is not under the base path
            }
        }
        
        if (relativePath.isEmpty()) {
            return false
        }
        
        // Count directory levels in the relative path
        int entryDepth = relativePath.split("/").length - 1
        
        // Include if within depth limit (depth 0 means only immediate children)
        return entryDepth <= depth
    }

    private RepositoryEntry createRepositoryEntry(TreeWalk treeWalk, String entryPath) {
        String name = entryPath.split('/').last()
        
        // Determine if it's a directory or file based on file mode
        EntryType type = treeWalk.isSubtree() ? EntryType.DIRECTORY : EntryType.FILE
        String sha = treeWalk.getObjectId(0).name()
        
        // For files, try to get size
        Long size = null
        if (type == EntryType.FILE) {
            try {
                def objectId = treeWalk.getObjectId(0)
                def loader = treeWalk.getObjectReader().open(objectId)
                size = loader.getSize()
            } catch (Exception e) {
                // Size not available, leave as null
            }
        }
        
        return new RepositoryEntry(
            name: name,
            path: entryPath,
            type: type,
            sha: sha,
            size: size
        )
    }

    @Override
    List<TagInfo> getTags() {
        final String prefix = 'refs/tags/'

        try ( final git = Git.open(new File(this.path, project))) {
            List<Ref> tags = git.tagList().call()
            final result = new ArrayList(tags.size())
            for (Ref ref : tags) {
                if( ref.name.startsWith(prefix) ) {
                    result.add( new TagInfo(ref.name.substring(prefix.length()), ref.getObjectId().name()))
                }
            }
            return result
        }
    }

    @Override
    List<BranchInfo> getBranches() {
        final String prefix = 'refs/heads/'

        try ( final git = Git.open(new File(this.path, project))) {
            List<Ref> tags = git.branchList().call()
            final result = new ArrayList(tags.size())
            for (Ref ref : tags) {
                if( ref.name.startsWith(prefix) ) {
                    result.add( new BranchInfo(ref.name.substring(prefix.length()), ref.getObjectId().name()) )
                }
            }
            return result
        }
    }

}
