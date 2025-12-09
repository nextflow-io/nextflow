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

package nextflow.scm

import groovy.transform.CompileStatic

/**
 * Context object that holds shared state and resources needed by repository strategies.
 * This avoids tight coupling between strategies and AssetManager by providing
 * a clean interface to access shared resources.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
class RepositoryContext {

    /**
     * The project name (e.g., "nextflow-io/hello")
     */
    final String project

    /**
     * Root directory where all projects are stored
     */
    final File root

    /**
     * The local root path for this project
     */
    final File localRootPath

    final RepositoryProvider provider


    RepositoryContext(
        String project,
        File root,
        File localRootPath,
        RepositoryProvider provider
    ) {
        this.project = project
        this.root = root
        this.localRootPath = localRootPath
        this.provider = provider
    }
}
