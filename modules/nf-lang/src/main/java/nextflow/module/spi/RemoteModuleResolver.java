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

import java.net.URI;
import java.nio.file.FileSystemNotFoundException;
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
     * Name of the directory under which modules are installed.
     */
    String MODULES_DIR = "modules";

    /**
     * Determine the base directory against which a remote module include should be resolved,
     * i.e. the directory whose {@code modules/} directory is searched.
     *
     * <p>An include in a vendored module's own script resolves against that module's directory,
     * so that a workflow module finds its dependencies in its own nested {@code modules/}
     * directory (nested vendoring). Any other script -- the entry script, or a plain local script
     * in the project -- resolves against the project directory, wherever in the project it lives.
     *
     * <p>A script is taken to belong to a vendored module when it lives under a {@code modules/}
     * directory within the project, which is where modules are installed at any depth.
     *
     * @param includingFile the script containing the include statement (may be null)
     * @param projectDir    the project directory (may be null when unknown, in which case the
     *                      including script's own directory is used)
     */
    static Path resolveBaseDir(Path includingFile, Path projectDir) {
        var parent = includingFile != null ? includingFile.getParent() : null;
        if( parent == null )
            return projectDir;
        if( projectDir == null )
            return parent;
        return isVendoredModuleDir(parent, projectDir) ? parent : projectDir;
    }

    /**
     * Variant of {@link #resolveBaseDir(Path, Path)} for a script identified by URI. A URI that
     * does not denote a file (e.g. an unsaved editor buffer in the language server) has no
     * directory of its own and resolves against the project directory.
     *
     * @param includingUri the URI of the script containing the include statement (may be null)
     * @param projectDir   the project directory
     */
    static Path resolveBaseDir(URI includingUri, Path projectDir) {
        if( includingUri == null )
            return projectDir;
        try {
            return resolveBaseDir(Path.of(includingUri), projectDir);
        }
        catch( IllegalArgumentException | FileSystemNotFoundException | SecurityException e ) {
            return projectDir;
        }
    }

    /**
     * @return true if the given directory is (or is nested within) a {@code modules/} directory
     * inside the project, i.e. it holds a vendored module rather than the project's own scripts.
     */
    private static boolean isVendoredModuleDir(Path directory, Path projectDir) {
        Path relative;
        try {
            relative = projectDir.toAbsolutePath().normalize()
                .relativize(directory.toAbsolutePath().normalize());
        }
        catch( IllegalArgumentException e ) {
            // not comparable to the project directory -- not a vendored module
            return false;
        }
        for( var segment : relative ) {
            var name = segment.toString();
            if( "..".equals(name) )
                return false;
            if( MODULES_DIR.equals(name) )
                return true;
        }
        return false;
    }

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
     * @param baseDir The directory relative to which the {@code modules/} directory is located:
     *                the project root for a top-level include, or the including module's own
     *                directory for a nested (vendored) include
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
