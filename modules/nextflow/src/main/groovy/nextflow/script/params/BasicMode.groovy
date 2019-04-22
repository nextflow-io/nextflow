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

import nextflow.script.TokenVar

/**
 * Process output mode
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
enum BasicMode implements OutParam.Mode {

    standard, flatten;

    static OutParam.Mode parseValue( def x ) {

        if( x instanceof OutParam.Mode ) {
            return x
        }

        def value = x instanceof TokenVar ? x.name : ( x instanceof String ? x : null )
        if( value ) {
            return BasicMode.valueOf(value)
        }

        throw new IllegalArgumentException("Not a valid output 'mode' value: $value")

    }
}