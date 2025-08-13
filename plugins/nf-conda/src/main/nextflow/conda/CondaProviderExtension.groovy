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

package nextflow.conda

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.packages.PackageProvider
import nextflow.packages.PackageProviderExtension

/**
 * Conda package provider extension
 * 
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
@CompileStatic
class CondaProviderExtension implements PackageProviderExtension {

    @Override
    PackageProvider createProvider(Session session) {
        def condaConfig = new CondaConfig(session.config.navigate('conda') as Map ?: [:])
        return new CondaPackageProvider(condaConfig)
    }

    @Override
    int getPriority() {
        return 100
    }
}