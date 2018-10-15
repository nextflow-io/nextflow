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
import groovy.util.logging.Slf4j
import static nextflow.Const.DEFAULT_BRANCH
import static nextflow.Const.DEFAULT_MAIN_FILE_NAME

/**
 * Models the nextflow config manifest settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class Manifest {

    private Map target

    Manifest() { target = Collections.emptyMap() }

    Manifest(Map object) {
        assert object != null
        this.target = new HashMap(object.size())
        final validFields = this.metaClass.properties.collect { it.name }.findAll { it!='class' }
        object.each { key, value ->
            if( validFields.contains(key) )
                target.put(key, value)
            else
                log.warn("Invalid config manifest attribute `$key`")
        }
    }

    String getHomePage() {
        target.homePage
    }


    String getDefaultBranch() {
        target.defaultBranch ?: DEFAULT_BRANCH
    }

    String getDescription() {
        target.description 
    }

    String getAuthor() {
        target.author
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

    String getName() {
        target.name
    }

}
