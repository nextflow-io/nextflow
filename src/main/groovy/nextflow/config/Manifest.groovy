/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

package nextflow.config

import groovy.transform.CompileStatic

import static nextflow.Const.DEFAULT_BRANCH
import static nextflow.Const.DEFAULT_MAIN_FILE_NAME

/**
 * Models the nextflow config manifest settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Manifest {

    private Map target

    Manifest() { target = Collections.emptyMap() }

    Manifest(Map object) {
        assert object != null
        this.target = object
    }

    String getHomePage() {
        target.homePage
    }


    String getDefaultBranch() {
        target.defaultBranch ?: DEFAULT_BRANCH
    }

    String getDescription() {
        // note: if description is not set it will return an empty ConfigObject
        // thus use the elvis operator to return null
        target.description ?: null
    }

    String getAuthor() {
        target.author ?: null
    }

    String getMainScript() {
        target.mainScript ?: DEFAULT_MAIN_FILE_NAME
    }

    def getGitmodules() {
        target.gitmodules
    }

    String getNextflowVersion() {
        target.nextflowVersion
    }

    String getVersion() {
        target.version
    }

}
