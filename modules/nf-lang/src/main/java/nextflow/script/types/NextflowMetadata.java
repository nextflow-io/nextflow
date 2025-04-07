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
package nextflow.script.types;

import nextflow.script.dsl.Constant;
import nextflow.script.dsl.Description;

public interface NextflowMetadata {

    @Constant("build")
    @Description("""
        Nextflow runtime build number.
    """)
    int getBuild();

    @Constant("timestamp")
    @Description("""
        Nextflow runtime compile timestamp.
    """)
    String getTimestamp();

    @Constant("version")
    @Description("""
        Nextflow runtime version number.
    """)
    VersionNumber getVersion();

}
