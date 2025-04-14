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
    A `MemoryUnit` represents a quantity of bytes.

    [Read more](https://nextflow.io/docs/latest/reference/stdlib.html#memoryunit)
""")
public interface MemoryUnit {

    @Description("""
        Get the memory value in bytes (B).
    """)
    long toBytes();

    @Description("""
        Get the memory value in gigabytes (rounded down).
    """)
    long toGiga();

    @Description("""
        Get the memory value in kilobytes (rounded down).
    """)
    long toKilo();

    @Description("""
        Get the memory value in megabytes (rounded down).
    """)
    long toMega();

    @Description("""
        Get the memory value in terms of a given unit (rounded down).
    """)
    long toUnit(String unit);

}
