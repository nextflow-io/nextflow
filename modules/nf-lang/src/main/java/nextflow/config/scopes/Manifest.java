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
package nextflow.config.scopes;

import java.util.List;
import java.util.Map;

import nextflow.config.schema.ConfigOption;
import nextflow.config.schema.ConfigScope;
import nextflow.script.dsl.Description;

public class Manifest implements ConfigScope {

    @ConfigOption
    @Description("""
        Project author name (use a comma to separate multiple names).
    """)
    public String author;

    @ConfigOption
    @Description("""
        List of project contributors. Should be a list of maps.
    """)
    public List<Map> contributors;

    @ConfigOption
    @Description("""
        Git repository default branch (default: `master`).
    """)
    public String defaultBranch;

    @ConfigOption
    @Description("""
        Free text describing the workflow project.
    """)
    public String description;

    @ConfigOption
    @Description("""
        Project documentation URL.
    """)
    public String docsUrl;

    @ConfigOption
    @Description("""
        Project related publication DOI identifier.
    """)
    public String doi;

    @ConfigOption
    @Description("""
        Project home page URL.
    """)
    public String homePage;

    @ConfigOption
    @Description("""
        Project related icon location (Relative path or URL).
    """)
    public String icon;

    @ConfigOption
    @Description("""
        Project license.
    """)
    public String license;

    @ConfigOption
    @Description("""
        Project main script (default: `main.nf`).
    """)
    public String mainScript;

    @ConfigOption
    @Description("""
        Project short name.
    """)
    public String name;

    @ConfigOption
    @Description("""
        Minimum required Nextflow version.
    """)
    public String nextflowVersion;

    @ConfigOption
    @Description("""
        Project organization.
    """)
    public String organization;

    @ConfigOption
    @Description("""
        Pull submodules recursively from the Git repository.
    """)
    public boolean recurseSubmodules;

    @ConfigOption
    @Description("""
        Project version number.
    """)
    public String version;

}
