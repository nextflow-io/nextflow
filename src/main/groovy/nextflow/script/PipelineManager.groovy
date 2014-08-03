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

package nextflow.script
import java.nio.file.Path

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.exception.AbortOperationException
import org.eclipse.jgit.api.CreateBranchCommand
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.api.errors.RefNotFoundException
/**
 * Handles operation on remote and local installed pipelines
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
class PipelineManager {

    static final String MANIFEST_FILE_NAME = '.MANIFEST'

    static final String SCRIPT_FILE_NAME = 'main.nf'

    static final String defaultOrganization = 'nextflow-io'

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

    private Git _git

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

    String resolveName( String name ) {
        assert name

        if( name.contains('/') ) {
            return name
        }

        def qualifiedName = find(name)
        if( !qualifiedName ) {
            return "$defaultOrganization/$name".toString()
        }

        if( qualifiedName instanceof List ) {
            throw new AbortOperationException("Which one do you mean?\n${qualifiedName.join('\n')}")
        }

        return qualifiedName
    }

    PipelineManager setPipeline( String name ) {
        assert name

        this.pipeline = resolveName(name)
        this.localPath = new File(root, pipeline)
        this.remoteUrl = new URL("https://api.github.com/repos/${pipeline}")
        this._git = null
        return this
    }

    PipelineManager setForce( boolean value ) {
        this.force = value
        return this
    }

    void checkValidGithubRepo() {

        def scriptName = scriptNameFor(remoteUrl)
        def scriptUrl = new URL("${remoteUrl}/contents/$scriptName")

        try {
            new JsonSlurper().parseText(scriptUrl.text)
        }
        catch( IOException e1 ) {

            try {
                new JsonSlurper().parseText(remoteUrl.text)
            }
            catch( IOException e2 ) {
                throw new AbortOperationException("Cannot find $pipeline pipeline -- Make sure exists a Github repository at http://github.com/$pipeline")
            }

            throw new AbortOperationException("Illegal pipeline repository http://github.com/$pipeline -- It must contain a script named '$SCRIPT_FILE_NAME' or a file '$MANIFEST_FILE_NAME' ")
        }

    }


    @Memoized
    String getGitRepositoryUrl() {

        if( localPath.exists() ) {
            return localPath.toURI().toString()
        }

        Map response
        try {
            response = new JsonSlurper().parseText(remoteUrl.text) as Map
        }
        catch( IOException e ) {
            throw new AbortOperationException("Cannot find a pipeline $pipeline -- Make sure exists a Github repository at http://github.com/$pipeline")
        }

        def result = response.get('clone_url')
        if( !result )
            throw new IllegalStateException("Missing clone URL for: $pipeline -- Github request: $remoteUrl")

        return result
    }

    File getLocalPath() { localPath }

    File getMainScriptFile() {
        if( !localPath.exists() ) {
            throw new AbortOperationException("Unknown pipeline folder: $localPath")
        }

        def mainScript = scriptNameFor(localPath)
        def result = new File(localPath, mainScript).absoluteFile
        if( !result.exists() )
            throw new AbortOperationException("Missing pipeline script: $result")

        return result
    }

    static String scriptNameFor( repository ) {
        def url
        if( repository instanceof URL )
            url = new URL("${repository}/contents/${MANIFEST_FILE_NAME}")
        else if( repository instanceof File || repository instanceof Path )
            url = new URL("file://${repository}/${MANIFEST_FILE_NAME}")
        else
            throw new IllegalArgumentException("Not a valid repository argument")

        def manifest = readManifest(url)
        return manifest?.get('main-script') ?: SCRIPT_FILE_NAME
    }


    static protected Map readManifest( URL url ) {
        InputStream input = null
        try {
            def manifest = new Properties()
            manifest.load( (input=url.openStream()) )
            return manifest
        }
        catch( IOException e ) {
            log.debug "Unable to read manifest file: $url"
            return null
        }
        finally {
            input?.close()
        }
    }

    String getBaseName() {
        pipeline.split('/')[1]
    }

    boolean isLocal() {
        localPath.exists()
    }

    /**
     * @return true if no differences exist between the working-tree, the index,
     *         and the current HEAD, false if differences do exist
     */
    boolean isClean() {
        git.status().call().isClean()
    }

    /**
     * @return The list of available pipelines
     */
    List<String> list() {
        log.debug "Listing pipelines in folders: $root"

        def result = []
        root.eachDir { File org ->
            org.eachDir { File it ->
                result << "${org.getName()}/${it.getName()}".toString()
            }
        }
        return result
    }

    def find( String name ) {
        def exact = []
        def partial = []

        list().each {
            def items = it.split('/')
            if( items[1] == name )
                exact << it
            else if( items[1].startsWith(name ) )
                partial << it
        }

        def list = exact ?: partial
        return list.size() ==1 ? list[0] : list
    }


    protected Git getGit() {
        if( !_git ) {
            _git = Git.open(localPath)
        }
        return _git
    }

    /**
     * Download a pipeline from a remote Github repository
     *
     * @param options
     */
    def download(String revision=null) {
        assert pipeline

        /*
         * if the pipeline already exists locally pull it from the remote repo
         */
        if( !localPath.exists() ) {
            // make sure it contains a valid repository
            checkValidGithubRepo()

            log.debug "Pulling $pipeline -- Using remote clone url: ${getGitRepositoryUrl()}"

            // clone it
            Git.cloneRepository()
                    .setURI(getGitRepositoryUrl())
                    .setDirectory(localPath)
                    .call()

        }


        log.debug "Pull pipeline $pipeline -- Using local path: $localPath"

        // verify that is clean
        if( !isClean() )
            throw new AbortOperationException("$pipeline contains uncommitted changes -- cannot pull from repository")

        if( revision && revision != getCurrentRevision() ) {
            /*
             * check out a revision before the pull operation
             */
            try {
                git.checkout() .setName(revision) .call()
            }
            /*
             * If the specified revision does not exist
             * Try to checkout it from a remote branch
             */
            catch ( RefNotFoundException e ) {
                return checkoutRemoteBranch(revision)
            }
        }

        // now pull to update it
        def result = git.pull().call()
        if(!result.isSuccessful())
            throw new AbortOperationException("Cannot pull pipeline: '$pipeline' -- ${result.toString()}")

        return result

    }

    /**
     * Clone a pipeline from a remote pipeline repository to the specified folder
     *
     * @param directory The folder when the pipeline will be cloned
     * @param revision The revision to be cloned. It can be a branch, tag, or git revision number
     */
    void clone(File directory, String revision = null) {

        def clone = Git.cloneRepository()
        def uri = getGitRepositoryUrl()
        log.debug "Clone pipeline $pipeline -- Using remote URI: ${uri} into: $directory"

        if( !uri )
            throw new AbortOperationException("Cannot find the specified pipeline: $pipeline")

        clone.setURI(uri)
        clone.setDirectory(directory)
        if( revision ) {
            clone.setBranch(revision)
        }

        clone.call()
    }

    String getCurrentRevision() {
        git.getRepository().getBranch()
    }

    void checkout( String revision = null ) {
        assert localPath

        def current = getCurrentRevision()
        if( current != masterRevision ) {
            if( !revision ) {
                throw new AbortOperationException("Pipeline '$pipeline' currently is sticked on revision: $current -- you need to specify explicitly a revision with the option '-r' to use it")
            }
        }
        else if( !revision || revision == current ) {
            // nothing to do
            return
        }

        // verify that is clean
        if( !isClean() )
            throw new AbortOperationException("Pipeline '$pipeline' contains uncommitted changes -- Cannot switch to revision: $revision")

        try {
            git.checkout().setName(revision) .call()
        }
        catch( RefNotFoundException e ) {
            checkoutRemoteBranch(revision)
        }

    }


    protected checkoutRemoteBranch( String revision ) {

        git.fetch().call()
        git.checkout()
                .setCreateBranch(true)
                .setName(revision)
                .setUpstreamMode(CreateBranchCommand.SetupUpstreamMode.TRACK)
                .setStartPoint("origin/" + revision)
                .call()
    }


}
