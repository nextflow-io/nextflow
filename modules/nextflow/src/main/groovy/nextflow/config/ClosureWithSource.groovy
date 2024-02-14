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

package nextflow.config

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope

/**
 * Placeholder class that contains a closure and its string
 * representation, which can be unwrapped to one or the other
 * based on runtime conditions.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@EqualsAndHashCode
@PackageScope
class ClosureWithSource {

    private Closure target

    private String str

    ClosureWithSource(Closure target, String str) {
        this.target = target
        this.str = str
    }

    Closure getTarget() { target }

    @Override
    String toString() { str }
}
