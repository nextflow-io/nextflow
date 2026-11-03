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
package nextflow.script.namespaces;

import groovy.lang.Closure;
import nextflow.script.dsl.Description;
import nextflow.script.dsl.Namespace;

public interface LogNamespace extends Namespace {

    @Description("""
        Log an error message to the console.
    """)
    void error(String message);

    @Description("""
        Log an info message to the console.
    """)
    void info(String message);

    @Description("""
        Log a warning message to the console.
    """)
    void warn(String message);

}
