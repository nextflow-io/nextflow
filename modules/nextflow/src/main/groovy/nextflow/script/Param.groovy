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

package nextflow.script

import java.lang.reflect.Type

import groovy.transform.Canonical
import groovy.transform.CompileStatic
/**
 * Models a declared parameter -- either a param declaration in the
 * `params` block or a `take:` input of a named workflow, which is
 * mapped from a pipeline parameter when the workflow is executed
 * directly.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Canonical
@CompileStatic
class Param {

    String name

    /** The declared type, or null if the param is untyped. */
    Type type

    /** Whether the declared type is nullable. */
    boolean optional

    /** The declared default value, or null. Workflow takes cannot declare a default. */
    Object defaultValue

}
