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
package nextflow.processor.hash

/**
 * Extracts the value(s) one hash key contributes for a task.
 *
 * A contributor chooses WHICH objects enter the hash; EncodingRules decide HOW any
 * object becomes bytes. Emitting the wrong NUMBER of values breaks byte-exactness
 * just as surely as emitting the wrong value, so an absent key must emit an empty
 * list — never a list containing null.
 */
interface Contributor {

    /**
     * Stable identity of this extraction, contributing to the spec fingerprint.
     * Two specs that extract differently for the same key must differ here.
     */
    String canonicalName()

    List<Object> emit(HashContext ctx)
}
