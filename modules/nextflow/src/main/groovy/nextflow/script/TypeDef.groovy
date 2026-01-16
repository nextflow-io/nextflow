/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.script

import groovy.transform.CompileStatic

/**
 * Models a type definition that can be included from a script
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class TypeDef extends ComponentDef {

    private Class target

    private String alias

    TypeDef(Class target) {
        this(target, target.getSimpleName())
    }

    TypeDef(Class target, String alias) {
        this.target = target
        this.alias = alias
    }

    @Override
    String getType() { 'type' }

    @Override
    String getName() { alias }

    @Override
    Object invoke_a(Object[] args) {
        throw new UnsupportedOperationException()
    }

    @Override
    TypeDef clone() {
        return (TypeDef)super.clone()
    }

    @Override
    TypeDef cloneWithName(String name) {
        def result = clone()
        result.@alias = name
        return result
    }
}
