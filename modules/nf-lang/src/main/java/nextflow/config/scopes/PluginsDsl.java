/*
 * Copyright 2024-2025, Seqera Labs
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
package nextflow.config.scopes;

import nextflow.script.dsl.Description;
import nextflow.script.dsl.DslScope;

public interface PluginsDsl extends DslScope {

    @Description("""
        Specify a plugin to be used by the pipeline. The plugin id can be a name (e.g. `nf-hello`) or a name with a version (e.g. `nf-hello@0.5.0`).
    """)
    public void id(String value);

}
