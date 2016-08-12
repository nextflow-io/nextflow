/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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

import groovy.transform.CompileStatic
import groovy.transform.ToString
import nextflow.scm.AssetManager

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
     * The
     */
    AssetManager.RevisionInfo revisionInfo

    /**
     * (local) project repository path
     */
    Path localPath

    /**
     * @return Directory where the main script file is stored
     */
    Path getParent() { main?.parent }

    /**
     * @return Main script file content as a text string
     */
    String getText() { main?.text }

    /**
     * @return Repository commitId
     */
    String getCommitId() { revisionInfo?.commitId }

    /**
     * @return Repository tag or branch
     */
    String getRevision() { revisionInfo?.revision }

    ScriptFile( File file ) {
        assert file
        main = file.toPath().complete()
        localPath = main.parent
    }

}
