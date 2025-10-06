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
            
            // Normalize path using base class helper
            def normalizedPath = normalizePath(path)
            
            if (normalizedPath && !normalizedPath.isEmpty()) {
                // Navigate to the specific directory first
                def dirWalk = TreeWalk.forPath(repo, normalizedPath, tree)
                try {
                    if (!dirWalk || !dirWalk.isSubtree()) {
                        return [] // Path doesn't exist or is not a directory
                    }
                    treeWalk.addTree(dirWalk.getObjectId(0))
                } finally {
                    dirWalk?.close()
                }
            } else {
                treeWalk.addTree(tree)
            }
            
            // For depth filtering, we need to traverse recursively when depth > 1
            // The shouldIncludeAtDepth filter will handle the actual depth limiting
            treeWalk.setRecursive(depth != 1)
            
            List<RepositoryEntry> entries = []
            
            while (treeWalk.next()) {
                String entryPath = treeWalk.getPathString()
                
                // Build full path for entries (relative paths need to be prefixed with base path)
                String fullPath = normalizedPath && !normalizedPath.isEmpty() ? "/" + normalizedPath + "/" + entryPath : "/" + entryPath
                
                // Filter by depth using base class helper
                if (shouldIncludeAtDepth(fullPath, path, depth)) {
                    entries.add(createRepositoryEntry(treeWalk, fullPath))
                }
            }
            
            treeWalk.close()
            revWalk.close()
            
            return entries.sort { it.name }
            
        } finally {
            git.close()
        }
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
