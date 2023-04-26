/*
 * Copyright 2013-2023, Seqera Labs
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
 * Models cpu resource request
 * 
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@CompileStatic
@EqualsAndHashCode
class CpuResource {

    final Integer request
    final Integer limit

    CpuResource( Number value ) {
        this(request: value)
    }

    CpuResource( Map res ) {
        if( res.request != null && res.limit != null ) {
            this.request = res.request as int
            this.limit = res.limit as int
        }
        else if( res.request != null ) {
            this.request = res.request as int
        }
        else if( res.limit != null ) {
            this.limit = res.limit as int
        }
    }

}
