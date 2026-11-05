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

package nextflow.cli

import groovy.transform.CompileStatic

/**
 * Web URLs of registry entries
 */
@CompileStatic
class RegistryUrls {

    /** Web URL of a registry entry; {@code kind} is the path segment, e.g. 'modules' or 'pipelines' */
    static String build(String registryUrl, String kind, String fullName, String version) {
        // Strip /api suffix to get the base UI URL
        final baseUrl = registryUrl.endsWith('/api') ? registryUrl[0..-5] : registryUrl
        return "${baseUrl}/${kind}/${fullName}@${version}"
    }
}
