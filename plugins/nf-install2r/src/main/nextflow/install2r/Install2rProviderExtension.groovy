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

package nextflow.install2r

import groovy.transform.CompileStatic
import nextflow.ISession
import nextflow.packages.PackageProvider
import nextflow.packages.PackageProviderExtension

/**
 * install2.r package provider extension
 */
@CompileStatic
class Install2rProviderExtension implements PackageProviderExtension {

    @Override
    PackageProvider createProvider(ISession session) {
        return new Install2rPackageProvider(new Install2rConfig(session.config.navigate('install2r') as Map ?: [:], System.getenv()))
    }

    @Override
    int getPriority() {
        return 100
    }
}
