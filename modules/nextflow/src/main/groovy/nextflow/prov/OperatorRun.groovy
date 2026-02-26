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
 *
 */

package nextflow.prov

import groovy.transform.Canonical
import groovy.transform.CompileStatic

/**
 * Model an operator run
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
@CompileStatic
class OperatorRun implements ProvLink {
    /**
     * The list of (object) ids that was received as input by a operator run
     */
    Set<Integer> inputIds = new LinkedHashSet<>(10)

    @Override
    String toString() {
        "OperatorRun[id=${System.identityHashCode(this)}; inputs=${inputIds}]"
    }
}
