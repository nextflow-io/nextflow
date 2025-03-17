/*
 * Copyright 2024-2025, Seqera Labs
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
package nextflow.script.types;

import java.util.List;
import java.util.Map;

import nextflow.script.dsl.Constant;
import nextflow.script.dsl.Description;

public interface Manifest {

    @Constant("author")
    @Description("""
        Project author name (use a comma to separate multiple names).
    """)
    String getAuthor();

    @Constant("contributors")
    @Description("""
        List of project contributors. Should be a list of maps.
    """)
    List<Map> getContributors();

    @Constant("defaultBranch")
    @Description("""
        Git repository default branch (default: `master`).
    """)
    String getDefaultBranch();

    @Constant("description")
    @Description("""
        Free text describing the workflow project.
    """)
    String getDescription();

    @Constant("docsUrl")
    @Description("""
        Project documentation URL.
    """)
    String getDocsUrl();

    @Constant("doi")
    @Description("""
        Project related publication DOI identifier.
    """)
    String getDoi();

    @Constant("homePage")
    @Description("""
        Project home page URL.
    """)
    String getHomePage();

    @Constant("icon")
    @Description("""
        Project related icon location (Relative path or URL).
    """)
    String getIcon();

    @Constant("license")
    @Description("""
        Project license.
    """)
    String getLicense();

    @Constant("mainScript")
    @Description("""
        Project main script (default: `main.nf`).
    """)
    String getMainScript();

    @Constant("name")
    @Description("""
        Project short name.
    """)
    String getName();

    @Constant("nextflowVersion")
    @Description("""
        Minimum required Nextflow version.
    """)
    String getNextflowVersion();

    @Constant("organization")
    @Description("""
        Project organization.
    """)
    String getOrganization();

    @Constant("recurseSubmodules")
    @Description("""
        Pull submodules recursively from the Git repository.
    """)
    boolean getRecurseSubmodules();

    @Constant("version")
    @Description("""
        Project version number.
    """)
    String getVersion();

}
