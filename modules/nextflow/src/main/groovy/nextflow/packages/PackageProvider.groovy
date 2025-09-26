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

import java.nio.file.Path

import groovy.transform.CompileStatic

/**
 * Interface for package providers (conda, pixi, mamba, etc.)
 * 
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
@CompileStatic
interface PackageProvider {

    /**
     * @return The name of this package provider (e.g., "conda", "pixi")
     */
    String getName()

    /**
     * @return Whether this package provider is available on the system
     */
    boolean isAvailable()

    /**
     * Create or get a cached environment for the given package specification
     * 
     * @param spec The package specification
     * @return The path to the environment
     */
    Path createEnvironment(PackageSpec spec)

    /**
     * Get the shell activation script for the given environment
     * 
     * @param envPath Path to the environment
     * @return Shell script snippet to activate the environment
     */
    String getActivationScript(Path envPath)

    /**
     * Check if the given package specification is valid for this provider
     * 
     * @param spec The package specification
     * @return True if the spec is valid for this provider
     */
    boolean supportsSpec(PackageSpec spec)

    /**
     * Get provider-specific configuration
     * 
     * @return Configuration object for this provider
     */
    Object getConfig()
}