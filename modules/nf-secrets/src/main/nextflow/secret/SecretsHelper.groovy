/*
 * Copyright 2021, Sage-Bionetworks
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

package nextflow.secret

import java.util.regex.Pattern

/**
 * Implements Secrets helper methods
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SecretsHelper {

    private final static Pattern NAME_REGEX = ~/^[a-zA-Z_](?:[0-9A-Za-z]+|(_)(?!\1)){1,49}$/

    static void checkName(String name) {
        final msg = "Invalid secret name: $name — It can only contains alphanumeric and non-consecutive underscore characters and cannot start with a numeric character"
        if( !NAME_REGEX.matcher(name).matches() )
            throw new IllegalArgumentException(msg)
        if( name.startsWith('__') )
            throw new IllegalArgumentException(msg)
    }

}
