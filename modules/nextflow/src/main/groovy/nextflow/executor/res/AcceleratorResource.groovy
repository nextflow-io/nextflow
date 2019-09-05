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

package nextflow.executor.res

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Models accelerator resource request
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@CompileStatic
@EqualsAndHashCode
class AcceleratorResource {

    final Integer request
    final Integer limit
    final String type
    final String runtime

    AcceleratorResource( Number value ) {
        this(limit: value)
    }

    AcceleratorResource( Map res ) {
        if( res.limit!=null && res.request!=null ) {
            this.limit = res.limit as int
            this.request = res.request as int
        }
        else if( res.limit!=null ) {
            this.limit = res.limit as int
            this.request = res.limit as int
        }
        else if( res.request != null ) {
            this.request = res.request as int
        }

        if( res.type )
            this.type = res.type as String

        if( res.runtime )
            this.runtime = res.runtime as String
    }

}
