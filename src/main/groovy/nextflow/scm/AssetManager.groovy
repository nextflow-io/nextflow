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

package nextflow.scm

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.exception.AbortOperationException
import org.eclipse.jgit.api.CreateBranchCommand
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.api.errors.RefNotFoundException
import org.eclipse.jgit.lib.Constants
import org.eclipse.jgit.lib.ObjectId
import org.eclipse.jgit.lib.Ref
import org.eclipse.jgit.lib.Repository
import org.eclipse.jgit.transport.UsernamePasswordCredentialsProvider

/**
 * Handles operation on remote and local installed pipelines
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
class AssetManager {

    static final String MANIFEST_FILE_NAME = 'nextflow.config'

    static final String DEFAULT_MAIN_FILE_NAME = 'main.nf'

    static final String DEFAULT_BRANCH = 'master'

    static final String DEFAULT_ORGANIZATION = System.getenv('NXF_ORG') ?: 'nextflow-io'

    static final String DEFAULT_HUB = System.getenv('NXF_HUB') ?: 'github'


    /**
     * The folder all pipelines scripts are installed
     */
    private File root = System.getenv('NXF_ASSETS') ? new File(System.getenv('NXF_ASSETS')) : Const.APP_HOME_DIR.resolve('assets').toFile()

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

    private Git _git

    private String _mainScript

    private RepositoryProvider provider

    String hub = DEFAULT_HUB

    String user

    String pwd

    /**
     * Create a new asset manager object with default parameters
     */
    AssetManager() {
    }

    /**
     * Create a new asset manager with the specified pipeline name
     *
     * @param pipeline The pipeline to be managed by this manager e.g. {@code nextflow-io/hello}
     */
    AssetManager(String pipeline) {
        setPipeline(pipeline)
    }

    /**
     * Create a new asset managed with the specified parameters. Accepted parameters are:
     * <li>hub
     * <li>user
     * <li>pas
     *
     * @param args
     */
    AssetManager(Map args) {
        this.hub = args?.hub ?: DEFAULT_HUB
        this.user = args?.user
        this.pwd = args?.pwd

        if( args.pipeline ) {
            setPipeline(args.pipeline as String)
        }

    }



    AssetManager setRoot( File root ) {
        assert root
        this.root = root
        return this
    }

    String resolveName( String name ) {
        assert name

        String[] parts = name.split('/')
        def last = parts[-1]
        if( last.endsWith('.nf') || last.endsWith('.nxf') ) {
            if( parts.size()==1 )
                throw new AbortOperationException("Not a valid pipeline name: $name")

            if( parts.size()==2 ) {
                _mainScript = last
                parts = [ parts.first() ]
            }
            else {
                _mainScript = parts[2..-1].join('/')
                parts = parts[0..1]
            }
        }

        if( parts.size() == 2 ) {
            return parts.join('/')
        }
        else if( parts.size()>2 ) {
            throw new AbortOperationException("Not a valid pipeline name: $name")
        }
        else {
            name = parts[0]
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

    AssetManager setPipeline( String name ) {
        assert name

        this.pipeline = resolveName(name)
        this.provider = createHubProviderFor(hub)
        this.localPath = new File(root, pipeline)
        this._git = null
        return this
    }

    String getPipeline() { pipeline }

    @PackageScope
    RepositoryProvider createHubProviderFor(String name) {
        if( name == 'github')
            return new GithubRepositoryProvider(pipeline: pipeline, user: user, pwd: pwd)

        if( name == 'bitbucket')
            return new BitbucketRepositoryProvider(pipeline: pipeline, user: user, pwd: pwd)

        throw new AbortOperationException("Unkwnon pipeline repository provider: $name")
    }

    AssetManager setLocalPath(File path) {
        this.localPath = path
        return this
    }

    AssetManager setForce( boolean value ) {
        this.force = value
        return this
    }

    void checkValidRemoteRepo() {
        def scriptName = getMainScriptName()
        provider.validateFor(scriptName)
    }


    @Memoized
    String getGitRepositoryUrl() {

        if( localPath.exists() ) {
            return localPath.toURI().toString()
        }

        provider.getCloneUrl()
    }

    File getLocalPath() { localPath }

    File getMainScriptFile() {
        if( !localPath.exists() ) {
            throw new AbortOperationException("Unknown pipeline folder: $localPath")
        }

        def mainScript = getMainScriptName()
        def result = new File(localPath, mainScript)
        if( !result.exists() )
            throw new AbortOperationException("Missing pipeline script: $result")

        return result
    }

    String getMainScriptName() {
        if( _mainScript )
            return _mainScript

        readManifest().mainScript ?: DEFAULT_MAIN_FILE_NAME
    }

    String getHomePage() {
        def manifest = readManifest()
        manifest.homePage ?: provider.getHomePage()
    }

    String getDefaultBranch() {
        readManifest().defaultBranch ?: DEFAULT_BRANCH
    }

    String getDescription() {
        // note: if description is not set it will return an empty ConfigObject
        // thus use the elvis operator to return null
        readManifest().description ?: null
    }

    protected Map readManifest() {
        ConfigObject result = null
        try {
            def text = localPath.exists() ? new File(localPath, MANIFEST_FILE_NAME).text : provider.readText(MANIFEST_FILE_NAME)
            if( text ) {
                def config = new ConfigSlurper().parse(text)
                result = (ConfigObject)config.manifest
            }
        }
        catch( Exception e ) {
            log.debug "Cannot read pipeline manifest -- Cause: ${e.message}"
        }
        // by default return an empty object
        return result ?: new ConfigObject()
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
            checkValidRemoteRepo()

            log.debug "Pulling $pipeline  -- Using remote clone url: ${getGitRepositoryUrl()}"

            // clone it
            def clone = Git.cloneRepository()
            if( user && pwd )
                clone.setCredentialsProvider( new UsernamePasswordCredentialsProvider(user, pwd) )

            clone
                .setURI(getGitRepositoryUrl())
                .setDirectory(localPath)
                .call()

            return "downloaded from ${gitRepositoryUrl}"
        }


        log.debug "Pull pipeline $pipeline  -- Using local path: $localPath"

        // verify that is clean
        if( !isClean() )
            throw new AbortOperationException("$pipeline  contains uncommitted changes -- cannot pull from repository")

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
        def pull = git.pull()
        if( user && pwd )
            pull.setCredentialsProvider( new UsernamePasswordCredentialsProvider(user, pwd))

        def result = pull.call()
        if(!result.isSuccessful())
            throw new AbortOperationException("Cannot pull pipeline: '$pipeline ' -- ${result.toString()}")

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
        log.debug "Clone pipeline $pipeline  -- Using remote URI: ${uri} into: $directory"

        if( !uri )
            throw new AbortOperationException("Cannot find the specified pipeline: $pipeline ")

        clone.setURI(uri)
        clone.setDirectory(directory)
        if( user && pwd )
            clone.setCredentialsProvider( new UsernamePasswordCredentialsProvider(user, pwd))

        if( revision )
            clone.setBranch(revision)

        clone.call()
    }

    /**
     * @return The symbolic name of the current revision i.e. the current checked out branch or tag
     */
    String getCurrentRevision() {
        Ref head = git.getRepository().getRef(Constants.HEAD);
        if( !head )
            return '(unknown)'

        if( head.isSymbolic() )
            return Repository.shortenRefName(head.getTarget().getName())

        if( !head.getObjectId() )
            return '(unknown)'

        // try to resolve the the current object it to a tag name
        Map<ObjectId, String> names = git.nameRev().addPrefix( "refs/tags/" ).add(head.objectId).call()
        names.get( head.objectId ) ?: head.objectId.name()
    }

    /**
     * @return A list of existing branches and tags names. For example
     * <pre>
     *     * master (default)
     *       patch-x
     *       v1.0 (t)
     *       v1.1 (t)
     * </pre>
     *
     * The star character on the left highlight the current revision, the string {@code (default)}
     *  ticks that it is the default working branch (usually the master branch), while the string {@code (t)}
     *  shows that the revision is a git tag (instead of a branch)
     */
    List<String> getRevisions() {

        def current = getCurrentRevision()
        def master = getDefaultBranch()

        List<String> branches = git.branchList()
                    .call()
                    .findAll { it.name.startsWith('refs/heads/') }
                    .collect { Repository.shortenRefName(it.name) }
                    .collect { (it == current ? '* ' : '  ') + it + ( master == it ? ' (default)' : '')  }

        List<String> tags = git.tagList()
                .call()
                .findAll  { it.name.startsWith('refs/tags/') }
                .collect { Repository.shortenRefName(it.name) }
                .collect { (it == current ? '* ' : '  ') + it + ' [t]'}

        def result = new ArrayList<String>(branches.size() + tags.size())
        result.addAll(branches)
        result.addAll(tags)
        return result
    }

    /**
     * Checkout a specific revision
     * @param revision The revision to be checked out
     */
    void checkout( String revision = null ) {
        assert localPath

        def current = getCurrentRevision()
        if( current != defaultBranch ) {
            if( !revision ) {
                throw new AbortOperationException("Pipeline '$pipeline ' currently is sticked on revision: $current -- you need to specify explicitly a revision with the option '-r' to use it")
            }
        }
        else if( !revision || revision == current ) {
            // nothing to do
            return
        }

        // verify that is clean
        if( !isClean() )
            throw new AbortOperationException("Pipeline '$pipeline ' contains uncommitted changes -- Cannot switch to revision: $revision")

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
