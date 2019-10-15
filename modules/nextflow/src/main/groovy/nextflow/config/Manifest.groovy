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

    /**
     * Controls whether repository sub-modules should be cloned along with the main one.
     *
     * @return
     *      Either a boolean value, a list object submodule names or a comma separated string
     *      of sub-module names
     */
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

    Map toMap() {
        final result = new HashMap(10)
        result.author = getAuthor()
        result.defaultBranch = getDefaultBranch()
        result.description = getDescription()
        result.homePage = homePage
        result.gitmodules = getGitmodules()
        result.mainScript = getMainScript()
        result.version = getVersion()
        result.nextflowVersion = getNextflowVersion()
        return result
    }
}
