/*
 * Copyright 2020, Microsoft Corp
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

package nextflow.cloud.azure.file

import java.util.regex.Matcher
import java.util.regex.Pattern

import groovy.transform.Canonical
import nextflow.extension.Bolts

/**
 * Model Azure bucket URL components
 */
@Canonical
class AzStorageContainerParser {

    private static final Pattern AZ_REGEX = ~/^azb:\/\/([^:]+:)(\/[^?]*)(?:\?account=(.+))?$/

    String container
    String path
    String account

    AzStorageContainerParser(String container, String path, String account) {
        this.container=container
        this.path=path
        this.account=account
    }

    AzStorageContainerParser(Matcher m) {
        container = m.group(1)
        path = m.group(2)
        account = m.group(3)
    }

    /**
     * @return The same as  {@link #path} but without the starting `/`
     */
    String getKey() {
        path ? Bolts.stripStart(path, '/') : path
    }

    static AzStorageContainerParser parse(String url) {
        final m = AZ_REGEX.matcher(url)
        if( m.matches() ) {
            return new AzStorageContainerParser(m)
        }
        return null
    }
}
