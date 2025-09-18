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

package nextflow.packages

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Unified package specification that supports multiple package managers
 * 
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
@CompileStatic
@ToString(includeNames = true)
@EqualsAndHashCode
class PackageSpec {

    /**
     * The package provider type (conda, pixi, mamba, micromamba, etc.)
     */
    String provider

    /**
     * List of package entries (e.g., ["samtools=1.17", "bcftools=1.18"])
     */
    List<String> entries

    /**
     * Environment file content (for YAML/TOML files)
     */
    String environment

    /**
     * Channels for package resolution (conda-specific)
     */
    List<String> channels

    /**
     * Additional provider-specific options
     */
    Map<String, Object> options

    PackageSpec() {
        this.entries = []
        this.channels = []
        this.options = [:]
    }

    PackageSpec(String provider, List<String> entries = [], Map<String, Object> options = [:]) {
        this.provider = provider
        this.entries = entries ?: []
        this.channels = []
        this.options = options ?: [:]
    }

    /**
     * Builder pattern methods
     */
    PackageSpec withProvider(String provider) {
        this.provider = provider
        return this
    }

    PackageSpec withEntries(List<String> entries) {
        this.entries = entries ?: []
        return this
    }

    PackageSpec withEnvironment(String environment) {
        this.environment = environment
        return this
    }

    PackageSpec withChannels(List<String> channels) {
        this.channels = channels ?: []
        return this
    }

    PackageSpec withOptions(Map<String, Object> options) {
        this.options = options ?: [:]
        return this
    }

    /**
     * Check if this spec is valid
     */
    boolean isValid() {
        return provider && (entries || environment)
    }

    /**
     * Check if this spec uses an environment file
     */
    boolean hasEnvironmentFile() {
        return environment != null && !environment.trim().isEmpty()
    }

    /**
     * Check if this spec has package entries
     */
    boolean hasEntries() {
        return entries && !entries.isEmpty()
    }
}