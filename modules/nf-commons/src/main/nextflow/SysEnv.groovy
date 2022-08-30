/*
 * Copyright 2020-2022, Seqera Labs
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

/**
 * Helper class that holds a reference system environment and
 * allow to swap to a different one for testing purposes
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class SysEnv {

    static final class Holder implements Map<String,String> {
        @Delegate Map<String,String> target
        Holder(Map<String,String> target) { this.target = target }
        void setTarget(Map<String,String> target) { this.target=target }
    }

    private static Holder holder = new Holder(System.getenv())

    static Map<String,String> get()  { return holder }

    static void set(Map<String,String> env) {
        holder.setTarget(env)
    }

    static void restore() {
        holder.setTarget(System.getenv())
    }

}
