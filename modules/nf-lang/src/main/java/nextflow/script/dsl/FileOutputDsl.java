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
package nextflow.script.dsl;

import java.nio.file.Path;
import java.util.Map;
import java.util.Set;

import groovy.transform.NamedParam;
import groovy.transform.NamedParams;

/**
 * The output functions that collect files from the task directory, shared by
 * {@link ProcessDsl.OutputDslV2} and {@link AgentDsl.AgentOutputDsl}.
 *
 * <p>Declared once because they mean the same thing in both: a process and an agent
 * lower {@code file(...)}/{@code files(...)} through the very same unstage visitor,
 * so an option added here must not be addable to only one of them.
 *
 * <p>Shared by the process and agent output scopes so the two cannot drift. Position in an
 * `extends` list does not matter: {@code VariableScopeChecker} searches every supertype.
 */
public interface FileOutputDsl extends DslScope {

    @Description("""
        Get a file from the task directory that matches the given pattern.
    """)
    Path file(
        @NamedParams({
            @NamedParam(value = "followLinks", type = Boolean.class),
            @NamedParam(value = "glob", type = Boolean.class),
            @NamedParam(value = "hidden", type = Boolean.class),
            @NamedParam(value = "includeInputs", type = Boolean.class),
            @NamedParam(value = "maxDepth", type = Integer.class),
            @NamedParam(value = "optional", type = Boolean.class),
            @NamedParam(value = "type", type = String.class),
        })
        Map<String,?> opts,
        String name
    );
    Path file(String name);

    @Description("""
        Get the files from the task directory that match the given pattern.
    """)
    Set<Path> files(
        @NamedParams({
            @NamedParam(value = "followLinks", type = Boolean.class),
            @NamedParam(value = "glob", type = Boolean.class),
            @NamedParam(value = "hidden", type = Boolean.class),
            @NamedParam(value = "includeInputs", type = Boolean.class),
            @NamedParam(value = "maxDepth", type = Integer.class),
            @NamedParam(value = "optional", type = Boolean.class),
            @NamedParam(value = "type", type = String.class),
        })
        Map<String,?> opts,
        String pattern
    );
    Set<Path> files(String pattern);

}
