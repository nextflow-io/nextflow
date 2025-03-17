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

@Description("""
    A `Duration` represents a duration of time.

    [Read more](https://nextflow.io/docs/latest/reference/stdlib.html#duration)
""")
public interface Duration {

    @Description("""
        Get the duration value in days (rounded down).
    """)
    long toDays();

    @Description("""
        Get the duration value in hours (rounded down).
    """)
    long toHours();

    @Description("""
        Get the duration value in milliseconds.
    """)
    long toMillis();

    @Description("""
        Get the duration value in minutes (rounded down).
    """)
    long toMinutes();

    @Description("""
        Get the duration value in seconds (rounded down).
    """)
    long toSeconds();

}
