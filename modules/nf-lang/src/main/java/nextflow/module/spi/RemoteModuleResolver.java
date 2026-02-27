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

package nextflow.module.spi;

import java.nio.file.Path;

/**
 * Service Provider Interface for resolving remote modules referenced with '@scope/name' syntax.
 *
 * <p>Implementations should handle:
 * <ul>
 *   <li>Checking if a module is already installed locally</li>
 *   <li>Downloading modules from a registry if not present</li>
 *   <li>Version resolution and validation</li>
 * </ul>
 *
 * <p>The interface follows the Java SPI pattern. Implementations should be registered
 * in META-INF/services/nextflow.module.spi.RemoteModuleResolver
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public interface RemoteModuleResolver {

    /**
     * Resolve a remote module reference (e.g., '@scope/name') to a local path.
     *
     * <p>This method should:
     * <ol>
     *   <li>Parse the module reference</li>
     *   <li>Check if the module is already installed locally</li>
     *   <li>Download and install the module if not present (auto-install)</li>
     *   <li>Validate version constraints if specified</li>
     * </ol>
     *
     * @param moduleName The module reference string (e.g., '@scope/name' or '@scope/name@version')
     * @param baseDir The base directory for the project (used to locate the modules directory)
     * @return Path to the resolved module's main.nf file
     * @throws IllegalArgumentException if the module reference is invalid or resolution fails
     */
    Path resolve(String moduleName, Path baseDir);

    /**
     * Get the priority of this resolver. Higher priority resolvers are tried first.
     *
     * <p>Use this to allow custom implementations to override the default resolver.
     * The default implementation should return 0. Custom implementations can return
     * positive values to take precedence.
     *
     * @return Priority value (higher = tried first), default should be 0
     */
    default int getPriority() {
        return 0;
    }
}