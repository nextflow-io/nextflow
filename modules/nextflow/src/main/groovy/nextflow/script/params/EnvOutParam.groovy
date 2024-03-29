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

package nextflow.script.params

import groovy.transform.InheritConstructors
import nextflow.script.TokenVar

/**
 * Model process `output: env PARAM` definition
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class EnvOutParam extends BaseOutParam implements OptionalParam {

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
        else if( obj instanceof CharSequence ) {
            this.nameObj = obj.toString()
        }
        else {
            throw new IllegalArgumentException("Unexpected environment output definition - it should be either a string or a variable identifier - offending value: ${obj?.getClass()?.getName()}")
        }

        return this
    }

}
