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

package nextflow.secret

import com.google.common.hash.Hasher
import groovy.transform.CompileStatic
import nextflow.util.CacheFunnel
import nextflow.util.CacheHelper

/**
 * Basic secret implementation
 *
 * It is declared as a record so that Gson can deserialize the secrets store file through the
 * canonical constructor. A regular class holding final fields would instead be created by the
 * Gson reflective adapter, which assigns the fields by reflection: as of Java 27 (JEP 500) that
 * emits a runtime warning and it is going to be rejected altogether in a future Java release.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
record SecretImpl(String name, String value) implements Secret, CacheFunnel {

    @Override
    String getName() { name }

    @Override
    String getValue() { value }

    @Override
    Hasher funnel(Hasher hasher, CacheHelper.HashMode mode) {
        hasher.putUnencodedChars(name)
        hasher.putUnencodedChars(value)
        return hasher
    }
}
