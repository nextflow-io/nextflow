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
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.exception.AbortOperationException
import nextflow.script.dsl.Description

import static nextflow.Const.DEFAULT_MAIN_FILE_NAME
/**
 * Models the nextflow config manifest settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("manifest")
@Description("""
    The `manifest` scope allows you to define some metadata that is useful when publishing or running your pipeline.
""")
@CompileStatic
class Manifest implements ConfigScope {

    @Deprecated
    @ConfigOption
    @Description("""
        Project author name (use a comma to separate multiple names).
    """)
    final String author

    @ConfigOption
    @Description("""
        List of project contributors. Should be a list of maps.
    """)
    final List<Contributor> contributors

    @ConfigOption
    @Description("""
        Git repository default branch (default: `master`).
    """)
    final String defaultBranch

    @ConfigOption
    @Description("""
        Free text describing the workflow project.
    """)
    final String description

    @ConfigOption
    @Description("""
        Project documentation URL.
    """)
    final String docsUrl

    @ConfigOption
    @Description("""
        Project related publication DOI identifier.
    """)
    final String doi

    @ConfigOption
    @Description("""
        Controls whether git sub-modules should be cloned with the main repository.

        Can be either a boolean value, a list of submodule names, or a comma-separated string of submodule names.
    """)
    final Object gitmodules

    @ConfigOption
    @Description("""
        Project home page URL.
    """)
    final String homePage

    @ConfigOption
    @Description("""
        Project related icon location (Relative path or URL).
    """)
    final String icon

    @ConfigOption
    @Description("""
        Project license.
    """)
    final String license

    @ConfigOption
    @Description("""
        Project main script (default: `main.nf`).
    """)
    final String mainScript

    @ConfigOption
    @Description("""
        Project short name.
    """)
    final String name

    @ConfigOption
    @Description("""
        Minimum required Nextflow version.
    """)
    final String nextflowVersion

    @ConfigOption
    @Description("""
        Project organization.
    """)
    final String organization

    @ConfigOption
    @Description("""
        Pull submodules recursively from the Git repository.
    """)
    final boolean recurseSubmodules

    @ConfigOption
    @Description("""
        Project version number.
    """)
    final String version

    /* required by extension point -- do not remove */
    Manifest() {}

    Manifest(Map opts) {
        author = opts.author as String
        contributors = parseContributors(opts.contributors)
        defaultBranch = opts.defaultBranch as String
        description = opts.description as String
        docsUrl = opts.docsUrl as String
        doi = opts.doi as String
        gitmodules = opts.gitmodules
        homePage = opts.homePage as String
        icon = opts.icon as String
        license = opts.license as String
        mainScript = opts.mainScript as String ?: DEFAULT_MAIN_FILE_NAME
        name = opts.name as String
        nextflowVersion = opts.nextflowVersion as String
        organization = opts.organization as String
        recurseSubmodules = opts.recurseSubmodules as boolean
        version = opts.version as String
    }

    private List<Contributor> parseContributors(Object value) {
        if( !value )
            return Collections.emptyList()

        try {
            final contributors = value as List<Map>
            return contributors.stream()
                .map(opts -> new Contributor(opts))
                .toList()
        }
        catch( IllegalArgumentException e ){
            throw new AbortOperationException(e.message)
        }
        catch( ClassCastException e ){
            throw new AbortOperationException("Invalid setting for `manifest.contributors` config option -- should be a list of maps")
        }
    }

    Map toMap() {
        return [
            author: author,
            contributors: contributors.stream().map(c -> c.toMap()).toList(),
            defaultBranch: defaultBranch,
            description: description,
            homePage: homePage,
            gitmodules: gitmodules,
            mainScript: mainScript,
            version: version,
            nextflowVersion: nextflowVersion,
            doi: doi,
            docsUrl: docsUrl,
            icon: icon,
            organization: organization,
            license: license,
        ]
    }

    @EqualsAndHashCode
    static class Contributor {
        final String name
        final String affiliation
        final String email
        final String github
        final List<ContributionType> contribution
        final String orcid

        Contributor(Map opts) {
            name = opts.name as String
            affiliation = opts.affiliation as String
            email = opts.email as String
            github = opts.github as String
            contribution = parseContributionTypes(opts.contribution as List<String>)
            orcid = opts.orcid as String
        }

        private List<ContributionType> parseContributionTypes(List<String> values) {
            if( values == null )
                return []
            final result = new LinkedList<ContributionType>()
            for( final value : values ) {
                try {
                    result.add(ContributionType.valueOf(value.toUpperCase()))
                }
                catch( IllegalArgumentException e ) {
                    throw new IllegalArgumentException("Invalid contribution type '$value' in `manifest.contributors` config option")
                }
            }
            return result.toSorted()
        }

        Map toMap() {
            return [
                name: name,
                affiliation: affiliation,
                email: email,
                github: github,
                contribution: contribution.stream().map(c -> c.toString().toLowerCase()).toList(),
                orcid: orcid,
            ]
        }
    }

    static enum ContributionType {
        AUTHOR,
        MAINTAINER,
        CONTRIBUTOR
    }
}
