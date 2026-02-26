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

import nextflow.config.schema.ConfigOption;
import nextflow.config.schema.ConfigScope;
import nextflow.script.dsl.Description;

public class TraceConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        When `true`, enables the creation of the execution trace report file (default: `false`).
    """)
    public boolean enabled;

    @ConfigOption
    @Description("""
        Comma separated list of fields to be included in the report.

        [Read more](https://nextflow.io/docs/latest/tracing.html#trace-report)
    """)
    public String fields;

    @ConfigOption
    @Description("""
        Trace file name (default: `'trace-<timestamp>.txt'`).
    """)
    public String file;

    @ConfigOption
    @Description("""
        When `true`, overwrites any existing trace file with the same name.
    """)
    public boolean overwrite;

    @ConfigOption
    @Description("""
        When `true`, uses raw number formatting i.e. durations are reported in milliseconds and memory in bytes.
    """)
    public boolean raw;

    @ConfigOption
    @Description("""
        Character used to separate values in each row (default: `\t`).
    """)
    public String sep;

}
