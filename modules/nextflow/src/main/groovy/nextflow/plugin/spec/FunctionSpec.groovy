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

import groovy.transform.CompileStatic
import nextflow.script.dsl.Description
import nextflow.script.types.Types
import org.codehaus.groovy.ast.ClassNode

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
                type: fromType(param.getType())
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

    private static Object fromType(Class c) {
        return fromType(new ClassNode(c))
    }

    private static Object fromType(ClassNode cn) {
        final name = Types.getName(cn.getTypeClass())
        if( !cn.isGenericsPlaceHolder() && cn.getGenericsTypes() != null ) {
            final typeArguments = cn.getGenericsTypes().collect { gt -> fromType(gt.getType()) }
            return [ name: name, typeArguments: typeArguments ]
        }
        else {
            return name
        }
    }
}
