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
package nextflow.script.types.shim;

import nextflow.script.dsl.Constant;
import nextflow.script.dsl.Description;
import nextflow.script.dsl.Ops;
import nextflow.script.types.Duration;
import nextflow.script.types.MemoryUnit;

@Description("""
    An integer is a whole number that can be positive or negative.

    [Read more](https://nextflow.io/docs/latest/reference/stdlib.html#integer)
""")
@Ops(IntegerOps.class)
public interface Integer {

    Integer intdiv(Integer b);

    // Duration

    @Constant("d")
    @Description("""
        Returns the equivalent Duration value in days.
    """)
    Duration d();

    @Constant("h")
    @Description("""
        Returns the equivalent Duration value in hours.
    """)
    Duration h();

    @Constant("m")
    @Description("""
        Returns the equivalent Duration value in minutes.
    """)
    Duration m();

    @Constant("ms")
    @Description("""
        Returns the equivalent Duration value in milliseconds.
    """)
    Duration ms();

    @Constant("s")
    @Description("""
        Returns the equivalent Duration value in seconds.
    """)
    Duration s();

    // MemoryUnit

    @Constant("B")
    @Description("""
        Returns the equivalent MemoryUnit value in bytes.
    """)
    MemoryUnit B();

    @Constant("KB")
    @Description("""
        Returns the equivalent MemoryUnit value in KB.
    """)
    MemoryUnit KB();

    @Constant("MB")
    @Description("""
        Returns the equivalent MemoryUnit value in MB.
    """)
    MemoryUnit MB();

    @Constant("GB")
    @Description("""
        Returns the equivalent MemoryUnit value in GB.
    """)
    MemoryUnit GB();

}

interface IntegerOps {

    Integer and(Integer a, Integer b);

    Integer bitwiseNegate();

    Integer compareTo(Integer a, Integer b);

    Integer compareTo(Integer a, Float b);

    Float div(Integer a, Integer b);

    Integer leftShift(Integer a, Integer b);

    Integer minus(Integer a, Integer b);

    Integer mod(Integer a, Integer b);

    Integer multiply(Integer a, Integer b);

    Integer negative();

    Integer or(Integer a, Integer b);

    Integer plus(Integer a, Integer b);

    Integer positive();

    Integer power(Integer a, Integer b);

    Integer rightShift(Integer a, Integer b);

    Integer xor(Integer a, Integer b);

}
