/*
 * Copyright 2020-2022, Seqera Labs
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
/**
 * Models a function component that can be included from a module script
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@CompileStatic
class FunctionDef extends ComponentDef implements ChainableDef {

    private BaseScript owner

    private String name

    private String alias

    FunctionDef(BaseScript owner, String name) {
        this.owner = owner
        this.name = name
        this.alias = name
    }

    protected FunctionDef() { }

    String getType() { 'function' }

    String getName() { alias }

    BaseScript getOwner() { owner }

    Object invoke_a(Object[] args) {
        final argsArr = ChannelOut.spread(args).toArray()
        final meta = owner.metaClass.getMetaMethod(name, argsArr)
        if( meta == null )
            throw new MissingMethodException(name, owner.getClass(), argsArr)
        Method callMethod = owner.getClass().getMethod(name, meta.getNativeParameterTypes())
        if( callMethod == null )
            throw new MissingMethodException(name, owner.getClass(), argsArr)
        return callMethod.invoke(owner, argsArr)
    }

    FunctionDef clone() {
        return (FunctionDef)super.clone()
    }

    FunctionDef cloneWithName(String name) {
        def result = clone()
        result.@alias = name
        return result
    }
}
