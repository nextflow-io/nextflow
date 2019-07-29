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

package nextflow.script

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.ToString
import nextflow.scm.AssetManager
import nextflow.util.CacheHelper

/**
 * Runnable pipeline script file
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true)
class ScriptFile {

    /**
     * Pipeline main script file
     */
    Path main

    /**
     * (remote) project repository URL
     */
    String repository

    /**
     * The revision information of the current version (commit id, branch/tag name)
     */
    AssetManager.RevisionInfo revisionInfo

    /**
     * (local) project repository path
     */
    Path localPath

    /**
     * The name of the project
     */
    String projectName 

    /**
     * @return Directory where the main script file is stored
     */
    Path getParent() { main?.parent }

    /**
     * @return Main script file content as a text string
     */
    String getText() { main?.text }

    /**
     * @return The repository commitId
     */
    String getCommitId() {
        revisionInfo?.commitId
    }

    String getScriptId() {
        main ? CacheHelper.hasher(main.text).hash().toString() : null
    }

    /**
     * @return Repository tag or branch
     */
    String getRevision() { revisionInfo?.name }

    ScriptFile( Path file ) {
        assert file
        main = file.complete()
        localPath = main.parent
    }

    ScriptFile( File file ) {
        this(file.toPath())
    }

}
