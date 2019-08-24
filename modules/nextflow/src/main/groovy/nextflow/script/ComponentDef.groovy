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

import groovy.transform.CompileStatic
import org.codehaus.groovy.runtime.InvokerHelper

/**
 * Models an abstract module component i.e. functions, processes
 * or (sub)workflow
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
abstract class ComponentDef implements Cloneable {

    abstract String getType()

    abstract String getName()

    abstract ComponentDef cloneWithName(String name)

    abstract Object invoke_a(Object[] args)

    Object invoke_o(Object args) {
        invoke_a(InvokerHelper.asArray(args))
    }

    String toString() {
        "${this.getClass().getSimpleName()}[$type $name]"
    }

    String getSignature() {
        "$type `$name`"
    }

}
