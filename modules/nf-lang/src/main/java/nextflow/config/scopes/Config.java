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
package nextflow.config.scopes;

import nextflow.config.spec.ConfigOption;
import nextflow.config.spec.ConfigScope;
import nextflow.script.dsl.Description;

public class Config implements ConfigScope {

    @Description("""
        The `env` scope allows you to define environment variables that will be exported into the environment where workflow tasks are executed.

        [Read more](https://nextflow.io/docs/latest/reference/config.html#env)
    """)
    public ConfigScope env;

    // NOTE: `nextflow` config options are inferred from FeatureFlagDsl
    public ConfigScope nextflow;

    @Description("""
        The `params` scope allows you to define parameters that will be accessible in the pipeline script.

        [Read more](https://nextflow.io/docs/latest/reference/config.html#params)
    """)
    public ConfigScope params;

    @ConfigOption
    @Description("""
        The `plugins` scope allows you to include plugins at runtime.

        [Read more](https://nextflow.io/docs/latest/plugins.html)
    """)
    public PluginsDsl plugins;

    // NOTE: `process` config options are inferred from ProcessDsl
    public ConfigScope process;

    @Description("""
        The `profiles` block allows you to define configuration profiles. A profile is a set of configuration settings that can be applied at runtime with the `-profile` command line option.

        [Read more](https://nextflow.io/docs/latest/config.html#config-profiles)
    """)
    public ConfigScope profiles;

}
