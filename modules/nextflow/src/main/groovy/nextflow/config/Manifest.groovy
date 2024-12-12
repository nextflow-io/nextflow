/*
 * Copyright 2013-2024, Seqera Labs
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

import java.util.stream.Collectors

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException

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
        target.defaultBranch
    }

    String getDefaultRevision() {
        target.defaultRevision ?: getVersion() ?: getDefaultBranch()
    }

    String getDescription() {
        target.description
    }

    String getAuthor() {
        target.author
    }

    List<Contributor> getContributors() {
        if( !target.contributors )
            return Collections.emptyList()

        try {
            final contributors = target.contributors as List<Map>
            return contributors.stream()
                .map(opts -> new Contributor(opts))
                .collect(Collectors.toList())
        }
        catch( ClassCastException | IllegalArgumentException e ){
            throw new AbortOperationException("Invalid config option `manifest.contributors` -- should be a list of maps")
        }
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

    boolean getRecurseSubmodules() {
        target.recurseSubmodules
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

    String getDoi() {
        target.doi
    }

    String getDocsUrl() {
        target.docsUrl
    }

    String getIcon() {
        target.icon
    }

    String getOrganization() {
        target.organization
    }

    String getLicense() {
        target.license
    }

    Map toMap() {
        final result = new HashMap(15)
        result.author = getAuthor()
        result.contributors = getContributors().stream()
            .map(c -> c.toMap())
            .collect(Collectors.toList())
        result.defaultBranch = getDefaultBranch()
        result.defaultRevision = getDefaultRevision()
        result.description = getDescription()
        result.homePage = homePage
        result.gitmodules = getGitmodules()
        result.mainScript = getMainScript()
        result.version = getVersion()
        result.nextflowVersion = getNextflowVersion()
        result.doi = getDoi()
        result.docsUrl = getDocsUrl()
        result.icon = getIcon()
        result.organization = getOrganization()
        result.license = getLicense()
        return result
    }

    @EqualsAndHashCode
    static class Contributor {
        String name
        String affiliation
        String email
        String github
        Set<ContributionType> contribution
        String orcid

        Contributor(Map opts) {
            name = opts.name as String
            affiliation = opts.affiliation as String
            email = opts.email as String
            github = opts.github as String
            contribution = (opts.contribution as List<String>).stream()
                .map(c -> ContributionType.valueOf(c.toUpperCase()))
                .collect(Collectors.toSet())
            orcid = opts.orcid as String
        }

        Map toMap() {
            final result = new HashMap(6)
            result.name = name
            result.affiliation = affiliation
            result.email = email
            result.github = github
            result.contribution = contribution.stream()
                .map(c -> c.toString().toLowerCase())
                .sorted()
                .collect(Collectors.toList())
            result.orcid = orcid
            return result
        }
    }

    static enum ContributionType {
        AUTHOR,
        MAINTAINER,
        CONTRIBUTOR
    }
}
