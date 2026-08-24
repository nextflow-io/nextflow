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
package nextflow.plugin.spec

import java.lang.reflect.Method
import java.lang.reflect.ParameterizedType
import java.lang.reflect.Type

import groovy.transform.CompileStatic
import nextflow.script.dsl.Description
import nextflow.script.dsl.Types

/**
 * Generate specs for functions, channel factories, and operators.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class FunctionSpec {

    static Map<String,?> of(Method method, String type) {
        final name = method.getName()
        final description = method.getAnnotation(Description)?.value()
        final returnType = fromType(method.getReturnType())
        final parameters = method.getParameters().collect { param ->
            [
                name: param.getName(),
                type: fromType(param.getParameterizedType())
            ]
        }

        return [
            type: type,
            spec: [
                name: name,
                description: description,
                returnType: returnType,
                parameters: parameters
            ]
        ]
    }

    private static Object fromType(Type type) {
        if( type instanceof ParameterizedType ) {
            final name = Types.getName(type.getRawType())
            final typeArguments = type.getActualTypeArguments().collect { t -> fromType(t) }
            return [ name: name, typeArguments: typeArguments ]
        }

        return Types.getName(type)
    }
}
