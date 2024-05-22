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
 * Models a function component that can be included from a module script
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@CompileStatic
class FunctionDef extends ComponentDef implements ChainableDef {

    private Object target

    private String name

    private String alias

    FunctionDef(Object target, String name) {
        this(target, name, name)
    }

    FunctionDef(Object target, String name, String alias) {
        this.target = target
        this.name = name
        this.alias = alias
    }

    protected FunctionDef() { }

    @Override
    String getType() { 'function' }

    @Override
    String getName() { alias }

    @Override
    Object invoke_a(Object[] args) {
        final arr = ChannelOut.spread(args).toArray()
        final meta = target.metaClass.getMetaMethod(name, arr)
        if( meta == null )
            throw new MissingMethodException(name, target.getClass(), arr)
        final types = meta.getNativeParameterTypes()
        final callMethod = target.getClass().getMethod(name, types)
        if( callMethod == null )
            throw new MissingMethodException(name, target.getClass(), arr)
        return callMethod.invoke(target, coerce0(arr,types))
    }

    static protected Object[] coerce0(Object[] arr, Class[] types) {
        assert arr.length==types.length
        // when the argument is a GString and type declared is a String
        // force the evaluation of the string to invoke the target function
        for( int i=0; i<arr.length; i++ ){
            if( arr[i] instanceof GString && types[i]==String.class )
                arr[i] = ((GString)arr[i]).toString()
        }
        return arr
    }

    @Override
    FunctionDef clone() {
        return (FunctionDef)super.clone()
    }

    @Override
    FunctionDef cloneWithName(String name) {
        def result = clone()
        result.@alias = name
        return result
    }
}
