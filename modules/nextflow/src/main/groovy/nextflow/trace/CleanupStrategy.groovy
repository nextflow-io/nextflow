/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.trace

import groovy.transform.CompileStatic

/**
 * Strategies to automatically cleanup task directories.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
enum CleanupStrategy {
    LAZY(1),
    EAGER(2),
    AGGRESSIVE(3)

    final int level

    CleanupStrategy(int level) {
        this.level = level
    }

    static boolean isValid(CharSequence name) {
        if( !name )
            return false
        try {
            valueOf(name.toString().toUpperCase())
            return true
        }
        catch( IllegalArgumentException e ) {
            return false
        }
    }
}
