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

import groovy.transform.CompileStatic

import javax.mail.MethodNotSupportedException
import java.lang.reflect.Method

import static ChannelOut.spreadToArray

/**
 * Models a function component that can be included from a module script
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@CompileStatic
class FunctionNameDef extends ComponentDef implements ChainableDef {

    private BaseScript owner

    private String name

    private String alias

    FunctionNameDef(BaseScript owner, String name) {
        this.owner = owner
        this.name = name
        this.alias = name
    }

    protected FunctionNameDef() { }

    String getType() { 'function' }

    String getName() { name }

    BaseScript getOwner() { owner }

    Object invoke_a(Object[] args) {
        throw new MethodNotSupportedException()
    }

    Object invoke_method(Object args){
        final Object[] argsArr = args instanceof Object[] ?  args as Object[]
                : args instanceof ChannelOut ? (args as List) as Object[]
                : new Object[]{args};
        def meta = owner.metaClass.getMetaMethod(name, argsArr)
        if( meta == null )
            throw new MissingMethodException(name, owner.getClass(), args)
        Method callMethod = owner.getClass().getMethod(name, meta.getNativeParameterTypes())
        if( callMethod == null )
            throw new MissingMethodException(name, owner.getClass(), args)
        callMethod.invoke(owner, argsArr)
    }

    FunctionNameDef clone() {
      return (FunctionNameDef)super.clone()
    }

    FunctionNameDef cloneWithName(String name) {
        def result = clone()
        result.@alias = name
        return result
    }
}
