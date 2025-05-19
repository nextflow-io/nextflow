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

import groovy.lang.Closure;
import nextflow.config.schema.ConfigOption;
import nextflow.config.schema.ConfigScope;
import nextflow.script.dsl.Description;

public class WorkflowConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        When `true`, the pipeline will exit with a non-zero exit code if any failed tasks are ignored using the `ignore` error strategy.
    """)
    public boolean failOnIgnore;

    @ConfigOption
    @Description("""
        Specify a closure that will be invoked at the end of a workflow run (including failed runs).
    """)
    public Closure onComplete;

    @ConfigOption
    @Description("""
        Specify a closure that will be invoked if a workflow run is terminated.
    """)
    public Closure onError;

    public WorkflowOutputConfig output;

    @Description("""
        The `workflow.output` scope provides options for publishing workflow outputs.
    
        [Read more](https://nextflow.io/docs/latest/reference/config.html#workflow)
    """)
    public WorkflowOutputConfig workflowOutput;

}
