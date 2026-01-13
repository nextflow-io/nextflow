/*
 * Copyright 2013-2026, Seqera Labs
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

import groovy.transform.*
import groovy.util.logging.Slf4j
import nextflow.config.Manifest
import nextflow.exception.AbortScmOperationException
import nextflow.script.ScriptFile
import nextflow.util.TestOnly
import org.eclipse.jgit.api.errors.RefNotFoundException
import org.eclipse.jgit.lib.Ref

/**
 * Handles operation on remote and local installed pipelines.
 * Uses the strategy pattern to support different ways to manage local installations.
 * Current available {@link RepositoryStrategy}:
 * - {@link LegacyRepositoryStrategy}: This is the traditional approach where each project gets a full git clone.
 * - {@Link MultiRevisionRepositoryStrategy}: This approach allows multiple revisions to coexist efficiently by sharing objects
 * through a bare repository and creating lightweight clones for each commit.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
class NextflowAssetManager extends AssetManager{

    private Manifest _manifest

    @TestOnly
    NextflowAssetManager(RepositoryProvider p = null){
        super(p)
    }


    NextflowAssetManager(String pipelineName, HubOptions options = null){
        super(pipelineName, ProviderConfigFactory.getDefault(), options)
    }

    NextflowAssetManager(String pipelineName, String revision, HubOptions options = null){
        super(pipelineName, ProviderConfigFactory.getDefault(), revision, options)
    }

    @Override
    RepositoryFactoryLoader getRepositoryFactory(){
        return NextflowRepositoryFactoryLoader.instance
    }

    private Manifest getManifest0() {
        String text = null
        try {
            text = localPath.exists() ? new File(localPath, ScmConst.MANIFEST_FILE_NAME).text : provider.readText(ScmConst.MANIFEST_FILE_NAME)
        }
        catch( FileNotFoundException e ) {
            log.debug "Project manifest does not exist: ${e.message}"
        }
        catch( Exception e ) {
            log.warn "Cannot read project manifest -- Cause: ${e.message ?: e}", e
        }
        return Manifest.parse(text)

    }

    Manifest getManifest() {
        if( _manifest == null ) {
            _manifest = getManifest0()
        }
        return _manifest
    }

    ScriptFile getScriptFile(String scriptName = null) {

        def result = new ScriptFile(getMainScriptFile(scriptName))
        result.revisionInfo = getCurrentRevisionAndName()
        result.repository = getRepositoryUrl()
        result.localPath = localPath.toPath()
        result.projectName = project

        return result
    }

    File getMainScriptFile(String scriptName=null) {
        if( !localPath.exists() ) {
            throw new AbortScmOperationException("Unknown project folder: $localPath")
        }

        def mainScript = scriptName ?: getMainScriptName()
        def result = new File(localPath, mainScript)
        if( !result.exists() )
            throw new AbortScmOperationException("Missing project main script: $result")

        return result
    }

    void checkValidRemoteRepo(String revision=null) {
        // Configure the git provider to use the required revision as source for all needed remote resources:
        // - config if present in repo (nextflow.config by default)
        // - main script (main.nf by default)
        provider.revision = revision
        final scriptName = getMainScriptName()
        provider.validateFor(scriptName)
    }

    String getMainScriptName() {
        return mainScript ?: manifest.mainScript
    }

    @Override
    String getDefaultBranch() {
        return manifest.getDefaultBranch()
            ?: super.getDefaultBranch()
    }

    String getHomePage() {
        getManifest().getHomePage() ?: provider.getRepositoryUrl()
    }

    /**
     * Download a pipeline from a remote Github repository
     *
     * @param revision The revision to download
     * @param deep Optional depth for shallow clones
     * @return A message representing the operation result
     */
    String download( String revision = null, Integer deep = null) {
        // If it is a new download check is a valid repository
        if( !localPath.exists() ) {
            checkValidRemoteRepo(revision)
        }
        super.download(revision, deep, manifest.recurseSubmodules)
    }

    /**
     * Clone a pipeline from a remote pipeline repository to the specified folder
     *
     * @param directory The folder when the pipeline will be cloned
     * @param revision The revision to be cloned. It can be a branch, tag, or git revision number
     */
    void clone(File directory,  String revision = null, Integer deep=null) {
        super.clone(directory,revision, deep, manifest.recurseSubmodules)
    }

    /**
     * Checkout a specific revision, and fetch remote if not locally available.
     * @param recurseSubmodules
     * @param revision The revision to be checked out
     */
    void checkout( String revision = null) {
        try {
            tryCheckout(revision)
        }
        catch( RefNotFoundException e ) {
            checkoutRemoteBranch(revision)
        }

    }

    protected Ref checkoutRemoteBranch(String revision) {
        super.checkoutRemoteBranch(revision, manifest.recurseSubmodules)
    }

    void updateModules() {
        super.updateModules(manifest.recurseSubmodules, manifest.gitmodules )
    }

}
