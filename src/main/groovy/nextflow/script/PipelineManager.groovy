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

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.exception.AbortOperationException
import org.eclipse.jgit.api.CreateBranchCommand
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.api.errors.RefNotFoundException
import org.eclipse.jgit.lib.Ref

/**
 * Handles operation on remote and local installed pipelines
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
class PipelineManager {

    static final String MANIFEST_FILE_NAME = '.PIPELINE-INF'

    static final String DEFAULT_MAIN_FILE_NAME = 'main.nf'

    static final String DEFAULT_MASTER_BRANCH = 'master'

    static final String DEFAULT_ORGANIZATION = System.getenv('NXF_ORG') ?: 'nextflow-io'

    /**
     * The folder all pipelines scripts are installed
     */
    private File root = new File(Const.APP_HOME_DIR, 'assets')

    /**
     * The pipeline name. It must be in the form {@code username/repo} where 'username'
     * is a valid user name or organisation account, while 'repo' is the repository name
     * containing the pipeline code
     */
    private String pipeline

    /**
     * Directory where the pipeline is cloned (i.e downloaded)
     */
    private File localPath

    /**
     * GitHub remote API entry point e.g. {@code https://api.github.com/repos/<pipeline>
     *
     * @link https://developer.github.com/v3/repos/
     *
     */
    private URL remoteUrl

    /**
     * url of github repository where the project is hosted
     */
    private String githubPage

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
            return "$DEFAULT_ORGANIZATION/$name".toString()
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
        this.githubPage = "http://github.com/$pipeline"
        this._git = null
        return this
    }

    PipelineManager setLocalPath(File path) {
        this.localPath = path
        return this
    }

    PipelineManager setForce( boolean value ) {
        this.force = value
        return this
    }

    void checkValidGithubRepo() {

        def scriptName = getMainScriptName()
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

            throw new AbortOperationException("Illegal pipeline repository $githubPage -- It must contain a script named '$DEFAULT_MAIN_FILE_NAME' or a file '$MANIFEST_FILE_NAME' ")
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

        def mainScript = getMainScriptName()
        def result = new File(localPath, mainScript).absoluteFile
        if( !result.exists() )
            throw new AbortOperationException("Missing pipeline script: $result")

        return result
    }

    String getMainScriptName() {
        readManifest()?.get('main-script') ?: DEFAULT_MAIN_FILE_NAME
    }

    String getHomePage() {
        readManifest()?.get('home-url') ?: githubPage
    }

    String getMasterBranch() {
        readManifest()?.get('master-branch') ?: DEFAULT_MASTER_BRANCH
    }

    protected Map readManifest() {
        try {
            if( localPath.exists() ) {
                def result = new Properties()
                result.load( new FileReader(new File(localPath, MANIFEST_FILE_NAME)) )
                return result
            }
            else {
                readManifestFrom("${remoteUrl}/contents/${MANIFEST_FILE_NAME}")
            }

        }
        catch( IOException e ) {
            log.debug "Cannot read pipeline manifest -- Cause: ${e.message}"
            return null
        }
    }

    /**
     * Given a Github content url: 1) get content object 2) decide from base64 3) read the content as {@link Properties} object
     *
     * @link https://developer.github.com/v3/repos/contents/#get-contents
     *
     * @param url
     * @return
     */
    @Memoized
    static protected Map readManifestFrom( String url ) {
        InputStream input = null
        try {
            Map response = (Map)new JsonSlurper().parseText(new URL(url).text)
            def bytes = response.get('content')?.toString()?.decodeBase64()

            def result = new Properties()
            result.load( new ByteArrayInputStream(bytes) )
            return result
        }
        catch( IOException e ) {
            log.debug "Unable to read manifest file: $url"
            return null
        }
        finally {
            input?.close()
        }
    }

    String getName() {
        return pipeline
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
        if( !root.exists() )
            return result

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
     * @param revision The revision to download
     * @result A message representing the operation result
     */
    def download(String revision=null) {
        assert pipeline

        /*
         * if the pipeline already exists locally pull it from the remote repo
         */
        if( !localPath.exists() ) {
            localPath.parentFile.mkdirs()
            // make sure it contains a valid repository
            checkValidGithubRepo()

            log.debug "Pulling $pipeline -- Using remote clone url: ${getGitRepositoryUrl()}"

            // clone it
            Git.cloneRepository()
                    .setURI(getGitRepositoryUrl())
                    .setDirectory(localPath)
                    .call()

            return "downloaded from ${gitRepositoryUrl}"
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
             * Try to checkout it from a remote branch and return
             */
            catch ( RefNotFoundException e ) {
                def ref = checkoutRemoteBranch(revision)
                return "checkout-out at ${ref.getObjectId()}"
            }
        }

        // now pull to update it
        def result = git.pull().call()
        if(!result.isSuccessful())
            throw new AbortOperationException("Cannot pull pipeline: '$pipeline' -- ${result.toString()}")

        return result?.mergeResult?.mergeStatus?.toString()

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
        if( current != masterBranch ) {
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


    protected Ref checkoutRemoteBranch( String revision ) {

        git.fetch().call()
        git.checkout()
                .setCreateBranch(true)
                .setName(revision)
                .setUpstreamMode(CreateBranchCommand.SetupUpstreamMode.TRACK)
                .setStartPoint("origin/" + revision)
                .call()
    }


}
