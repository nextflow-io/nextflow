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
 *
 */

package nextflow

import groovy.transform.CompileStatic
import groovy.transform.PackageScope

/**
 * Helper class that holds a reference system environment and
 * allow to swap to a different one for testing purposes
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class SysEnv {

    private static Holder holder = new Holder(System.getenv())

    private static final List<Map<String,String>> history = new ArrayList<Map<String,String>>()

    static boolean containsKey(String key) {
        return holder.containsKey(key)
    }
    
    static Map<String,String> get()  {
        return holder
    }

    static String get(String name) {
        return holder.get(name)
    }

    static String get(String name, String defValue) {
        return holder.containsKey(name) ? holder.get(name) : defValue
    }

    static boolean getBool(String name, boolean defValue) {
        final result = get(name,String.valueOf(defValue))
        return Boolean.parseBoolean(result)
    }

    static void push(Map<String,String> env) {
        history.push(holder.getTarget())
        holder.setTarget(env)
    }

    static void pop() {
        holder.setTarget(history.pop())
    }

}

@PackageScope
class Holder implements Map<String,String> {
    @Delegate Map<String,String> target
    Holder(Map<String,String> target) { this.target = target }
    void setTarget(Map<String,String> target) { this.target=target }
}

