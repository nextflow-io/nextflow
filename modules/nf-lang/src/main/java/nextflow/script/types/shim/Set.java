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
    A set is an unordered collection that cannot contain duplicate elements.

    [Read more](https://nextflow.io/docs/latest/reference/stdlib.html#sete)
""")
@Ops(SetOps.class)
public interface Set<E> extends Iterable<E> {

    @Description("""
        Returns the intersection of the set and the given iterable.
    """)
    Set<E> intersect(Iterable<E> right);
}

interface SetOps<E> {

    Boolean isCase(E a, Set<E> b);

    Set<E> minus(Set<E> a, Iterable<E> b);

    Set<E> plus(Set<E> a, Iterable<E> b);

}
