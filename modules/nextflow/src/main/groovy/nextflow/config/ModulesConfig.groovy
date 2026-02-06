/*
 * Copyright 2013-2026, Seqera Labs
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
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.TestOnly

/**
 * Configuration scope for module version declarations
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@ScopeName("modules")
@Description("""
    The `modules` scope provides module version declarations for the Nextflow module system.
    Each entry maps a module reference to a specific version.
""")
@CompileStatic
class ModulesConfig implements ConfigScope {

    @ConfigOption
    @Description("Module version mappings (module name -> version)")
    private Map<String, String> modules = [:]

    /* required by extension point -- do not remove */
    ModulesConfig() {}

    ModulesConfig(Map opts) {
        if (opts) {
            opts.each { key, value ->
                this.modules.put(key.toString(), value.toString())
            }
        }
    }

    /**
     * Get the configured version for a module
     *
     * @param moduleName The module name (e.g., "@nf-core/fastqc")
     * @return The configured version, or null if not configured
     */
    String getVersion(String moduleName) {
        return this.modules.get(moduleName)
    }

    /**
     * Get all configured modules
     *
     * @return Map of module name to version
     */
    @TestOnly
    Map<String, String> getAllModules() {
        return Collections.unmodifiableMap(modules)
    }

    /**
     * Set a module version
     *
     * @param moduleName The module name
     * @param version The version to set
     */
    void setVersion(String moduleName, String version) {
        this.modules.put(moduleName, version)
    }

    /**
     * Check if a module version is configured
     *
     * @param moduleName The module name
     * @return true if version is configured
     */
    boolean hasVersion(String moduleName) {
        return modules.containsKey(moduleName)
    }
}
