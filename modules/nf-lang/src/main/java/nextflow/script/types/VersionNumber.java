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

import nextflow.script.dsl.Description;

public interface VersionNumber {

    @Description("""
        Get the major version number, i.e. the first version component.
    """)
    String getMajor();

    @Description("""
        Get the minor version number, i.e. the second version component.
    """)
    String getMinor();

    @Description("""
        Get the minor version number, i.e. the third version component.
    """)
    String getPatch();

    @Description("""
        Check whether the version satisfies a version requirement.
    """)
    boolean matches(String condition);

}
