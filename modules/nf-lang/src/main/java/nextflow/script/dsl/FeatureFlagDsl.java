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
package nextflow.script.dsl;

public class FeatureFlagDsl {

    @Deprecated
    @FeatureFlag("nextflow.enable.configProcessNamesValidation")
    @Description("""
        When `true`, prints a warning for every `withName:` process selector that doesn't match a process in the pipeline (default: `true`).
    """)
    public boolean configProcessNamesValidation;

    @Deprecated
    @FeatureFlag("nextflow.enable.dsl")
    @Description("""
        Defines the DSL version (`1` or `2`).
    """)
    public float dsl;

    @FeatureFlag("nextflow.enable.moduleBinaries")
    @Description("""
        When `true`, enables the use of modules with executable scripts i.e. [module binaries](https://nextflow.io/docs/latest/module.html#module-binaries).
    """)
    public boolean moduleBinaries;

    @FeatureFlag("nextflow.enable.strict")
    @Description("""
        When `true`, the pipeline is executed in [strict mode](https://nextflow.io/docs/latest/reference/feature-flags.html).
    """)
    public boolean strict;

    @FeatureFlag("nextflow.preview.output")
    @Description("""
        When `true`, enables the use of the [workflow output definition](https://nextflow.io/docs/latest/workflow.html#workflow-output-def).
    """)
    public boolean previewOutput;

    @FeatureFlag("nextflow.preview.recursion")
    @Description("""
        When `true`, enables the use of [process and workflow recursion](https://github.com/nextflow-io/nextflow/discussions/2521).
    """)
    public boolean previewRecursion;

    @FeatureFlag("nextflow.preview.package")
    @Description("""
        When `true`, enables the unified package management system with the `package` directive.
    """)
    public boolean previewPackage;

}
