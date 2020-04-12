/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import java.lang.reflect.Method

import groovy.transform.CompileStatic
import static ChannelOut.spreadToArray

/**
 * Models a function component that can be included from a module script
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@CompileStatic
class FunctionDef extends ComponentDef implements ChainableDef {

    private BaseScript owner

    private Method method

    private String name

    FunctionDef(BaseScript owner, Method method) {
        this.owner = owner
        this.method = method
        this.name = method.name
    }

    protected FunctionDef() { }

    String getType() { 'function' }

    Method getMethod() { method }

    String getName() { name }

    BaseScript getOwner() { owner }

    Object invoke_a(Object[] args) {
        method.invoke(owner, spreadToArray(args))
    }

    FunctionDef clone() {
      return (FunctionDef)super.clone()
    }

    FunctionDef cloneWithName(String name) {
        def result = clone()
        result.@name = name
        return result
    }
}
