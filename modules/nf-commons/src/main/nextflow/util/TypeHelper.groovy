/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.util

import java.lang.reflect.ParameterizedType
import java.lang.reflect.Type

import groovy.transform.CompileStatic

/**
 * A utility class that provides helper methods for working with generic types at runtime.
 * <p>
 * This class is designed to extract type information from objects that have generic superclasses.
 * </p>
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TypeHelper {

    /**
     * Retrieves the generic type at the specified index from the given object's superclass.
     *
     * <p>This method assumes that the object's class extends a parameterized type,
     * and it returns the type argument at the given index.</p>
     *
     * @param object the object whose generic type is to be retrieved
     * @param index  the index of the generic type parameter (starting from 0)
     * @return the {@link Type} representing the generic type at the specified index
     *
     * @example
     * <pre>
     * class ExampleClass extends GenericBase&lt;String, Integer&gt; {}
     *
     * Type type = TypeHelper.getGenericType(new ExampleClass(), 1);
     * System.out.println(type); // Output: class java.lang.Integer
     * </pre>
     */
    static Type getGenericType(Object object, int index) {
        final params = (ParameterizedType) (object.getClass().getGenericSuperclass())
        return params.getActualTypeArguments()[index]
    }

}
