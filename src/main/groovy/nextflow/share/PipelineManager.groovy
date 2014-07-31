/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.share
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.exception.InvalidRepositoryStateException
import org.eclipse.jgit.api.Git

/**
 * Handles operation on remote and local installed pipelines
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
class PipelineManager {

    static final String MAIN_FILE_NAME = 'main.nf'

    /**
     * The folder all pipelines scripts are installed
     */
    private File root = new File(Const.APP_HOME_DIR, 'repos')

    /**
     * The current pipeline name. It must be in the form {@code username/repo} where 'username'
     * is a valid user name or organisation account, while 'repo' is the repository name
     * containing the
     */
    private String pipeline

    private File localPath

    private URL remoteUrl

    private String masterRevision = 'master'

    PipelineManager() {
    }

    PipelineManager(String name) {
        setPipeline(name)
    }

    PipelineManager setRoot( File root ) {
        assert root
        this.root = root
        return this
    }

    PipelineManager setPipeline( String name ) {
        assert name
        this.pipeline = name
        this.localPath = new File(root, name)
        this.remoteUrl = new URL("https://api.github.com/repos/${name}")
        return this
    }

    String getRepoUri() {

        if( localPath.exists() ) {
            return localPath.toURI().toString()
        }

        def repo = new JsonSlurper().parseText(remoteUrl.text) as Map
        return repo.get('clone_url')
    }

    File getLocalPath() { localPath }

    String getRemoteUri() { remoteUrl.toURI().toString() }

    File getMainFile() {
        return new File(localPath, MAIN_FILE_NAME).absoluteFile
    }

    String getBaseName() {
        pipeline.split('/')[1]
    }

    boolean isLocal() {
        localPath.exists()
    }

    /**
     * @return The list of available pipelines
     */
    List<String> list() {
        log.debug "Listing pipelines in folders: $root"

        def result = []
        root.eachDir { File org ->
            org.eachDir { File it ->
                result << "${org.getName()}/${it.getName()}"
            }
        }
        return result
    }

    /**
     * Pull a pipeline from a remoteGithub repository
     *
     * @param options
     */
    def void pull() {
        assert pipeline

        /*
         * if the pipeline already exists locally pull it from the remote repo
         */
        if( localPath.exists() ) {
            log.debug "Pull pipeline $pipeline -- Using local path: $localPath"
            def git = Git.open(localPath)
            def pull = git.pull()

            def result = pull.call()
            if(!result.isSuccessful())
                throw new InvalidRepositoryStateException("Cannot pull pipeline: '$pipeline' -- ${result.toString()}")
        }
        /*
         * otherwise try to clone it from the remote repository
         */
        else {
            log.debug "Pull pipeline $pipeline -- Using remote URI: ${getRepoUri()}"
            def clone = Git.cloneRepository()
            clone.setURI(getRepoUri())
            clone.setDirectory( localPath )
            clone.call()
        }
    }

    /**
     * Clone a pipeline from a remote pipeline repository to a specified folder
     *
     * @param options
     */
    void clone(File directory, String revision = null) {

        def clone = Git.cloneRepository()
        def uri = getRepoUri()
        log.debug "Clone pipeline $pipeline -- Using remote URI: ${uri} into: $directory"

        if( !uri )
            throw new IllegalArgumentException("Cannot find the specified pipeline: $pipeline")

        clone.setURI(uri)
        clone.setDirectory(directory)

        if( revision ) {
            clone.setBranch(revision)
        }

        clone.call()
    }

    String getCurrentRevision() {
        def git = Git.open(localPath)
        try {
            return git.getRepository().getBranch()
        }
        finally {
            git.close()
        }
    }


    void checkout( String revision ) {
        assert localPath

        if( !revision )
            revision = masterRevision

        def git = Git.open(localPath)
        def current = git.getRepository().getBranch()
        if( current == revision ) {
            // nothing to do
            return
        }

        // verify that is clean
        def isClean =  git.status().call().isClean()
        if( !isClean )
            throw new IllegalStateException("Pipeline '$pipeline' contains uncommitted changes -- Cannot switch to revision: $revision")


        git.checkout().setName(revision) .call()

    }

    /**
     * @return true if no differences exist between the working-tree, the index,
     *         and the current HEAD, false if differences do exist
     */
    boolean isClean() {
        Git.open(localPath).status().call().isClean()
    }

}
