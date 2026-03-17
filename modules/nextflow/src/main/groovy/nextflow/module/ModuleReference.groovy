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

package nextflow.module

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import nextflow.exception.AbortOperationException

import java.util.regex.Pattern

/**
 * Represents a reference to a module in DSL include statements
 *
 * @author Jorge Ejarqe <jorge.ejarque@seqera.io>
 */
@CompileStatic
@EqualsAndHashCode
class ModuleReference {

    // Must stay in sync with ModuleResolver.REMOTE_MODULE_PATTERN in nf-lang
    private static final Pattern MODULE_NAME_PATTERN = ~/^([a-z0-9][a-z0-9._\-]*)\/([a-z][a-z0-9._\-]*(?:\/[a-z][a-z0-9._\-]*)*)$/

    final String scope
    final String name
    final String fullName

    ModuleReference(String scope, String name) {
        this.scope = scope
        this.name = name
        this.fullName = "${scope}/${name}"
    }

    /**
     * Parse a module reference from a string as "scope/name"
     *
     * @param source The module reference string
     * @return A ModuleReference object
     * @throws AbortOperationException if the format is invalid
     */
    static ModuleReference parse(String source) {
        if( !source ) {
            throw new AbortOperationException("Module reference cannot be empty")
        }

        // Trim whitespace
        source = source.trim()

        def matcher = MODULE_NAME_PATTERN.matcher(source)
        if( !matcher.matches() ) {
            throw new AbortOperationException(
                "Invalid module reference: '${source}'. " +
                    "Expected format: scope/name where scope is lowercase alphanumeric with dots/underscores/hyphens " +
                    "and name is lowercase alphanumeric with underscores/hyphens, optionally with slash-separated segments"
            )
        }

        return new ModuleReference(matcher.group(1), matcher.group(2))
    }

    @Override
    String toString() {
        return fullName
    }
}
