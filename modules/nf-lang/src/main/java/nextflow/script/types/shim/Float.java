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

import nextflow.script.dsl.Description;
import nextflow.script.dsl.Ops;

@Description("""
    A float is a floating-point number (i.e. real number) that can be positive or negative.

    [Read more](https://nextflow.io/docs/latest/reference/stdlib.html#float)
""")
@Ops(FloatOps.class)
public interface Float {
    Integer intValue();
}

interface FloatOps {

    Float and(Float a, Float b);

    Integer compareTo(Float a, Float b);

    Integer compareTo(Float a, Integer b);

    Float div(Float a, Float b);

    Float minus(Float a, Float b);

    Float multiply(Float a, Float b);

    Float negative();

    Float ofType(Integer n);

    Float plus(Float a, Float b);

    Float positive();

    Float power(Float a, Float b);

}
