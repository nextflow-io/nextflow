/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
    protected byte[] readBytes(String path) {

        final git = Git.open(new File(this.path, project))
        try {
            final repo = git.getRepository()

            def lastCommitId = repo.resolve(Constants.HEAD)
            def revWalk = new RevWalk(repo)
            def commit = revWalk.parseCommit(lastCommitId)
            def tree = commit.getTree()

            def treeWalk = TreeWalk.forPath(repo, path, tree)
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


}
