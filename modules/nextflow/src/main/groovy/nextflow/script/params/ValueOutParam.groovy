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

package nextflow.script.params

import groovy.transform.InheritConstructors
import nextflow.script.TokenVar


/**
 * Model a process *value* output parameter
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class ValueOutParam extends BaseOutParam {

    protected target

    String getName() {
        return nameObj ? super.getName() : null
    }


    BaseOutParam bind( def obj ) {
        // the target value object
        target = obj

        // retrieve the variable name to be used to fetch the value
        if( obj instanceof TokenVar ) {
            this.nameObj = obj.name
        }

        return this
    }

    /**
     * Given the {@link nextflow.processor.TaskContext} object resolve the actual value
     * to which this param is bound
     *
     * @param context An instance of {@link nextflow.processor.TaskContext} holding the task evaluation context
     * @return The actual value to which this out param is bound
     */
    def resolve( Map context ) {

        switch( target ) {
            case TokenVar:
                return context.get(target.name)

            case Closure:
                return target.cloneWith(context).call()

            case GString:
                return target.cloneAsLazy(context).toString()

            default:
                return target
        }
    }

}
