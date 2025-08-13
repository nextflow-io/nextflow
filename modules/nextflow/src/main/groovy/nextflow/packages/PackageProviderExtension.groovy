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
import nextflow.Session
import org.pf4j.ExtensionPoint

/**
 * Extension point for package provider plugins
 * 
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
@CompileStatic
interface PackageProviderExtension extends ExtensionPoint {

    /**
     * Create a package provider instance
     * 
     * @param session The current Nextflow session
     * @return A package provider instance
     */
    PackageProvider createProvider(Session session)

    /**
     * Get the priority of this extension (higher values take precedence)
     * 
     * @return Priority value
     */
    int getPriority()
}