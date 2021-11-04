/*
 * Copyright 2020-2021, Seqera Labs
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
 *
 */

package nextflow.sraql.config

import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Model a SRAQL dataSource configuration
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
@ToString(includePackage = false, includeNames = true)
@EqualsAndHashCode
class SraqlDataSource {
    public static String DEFAULT_SOURCE = 'google-bigquery'

    static SraqlDataSource DEFAULT = new SraqlDataSource(Collections.emptyMap())

    String source

    SraqlDataSource(Map opts) {
        this.source = opts.source ?: DEFAULT_SOURCE
    }

    Map toMap() {
        final result = new HashMap(10)

        if( source )
            result.source = source

        return result
    }

    @Override
    String toString() {
        return "SraqlDataSource[source=$source]"
    }
}
