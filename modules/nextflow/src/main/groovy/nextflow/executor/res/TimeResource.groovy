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

package nextflow.executor.res

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import nextflow.util.Duration

/**
 * Models time resource request
 * 
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@CompileStatic
@EqualsAndHashCode
class TimeResource {

    final Duration request
    final Duration limit

    TimeResource( value ) {
        this(limit: value)
    }

    TimeResource( Map res ) {
        if( res.request != null && res.limit != null ) {
            this.request = toDuration(res.request)
            this.limit = toDuration(res.limit)
        }
        else if( res.request != null ) {
            this.request = toDuration(res.request)
        }
        else if( res.limit != null ) {
            this.request = toDuration(res.limit)
            this.limit = toDuration(res.limit)
        }
    }

    Duration toDuration( value ) {
        if ( value instanceof Duration ) {
            return value as Duration
        }

        try {
            return new Duration(value.toString().trim())
        }
        catch( Exception e ) {
            throw new IllegalArgumentException("Not a valid duration value: $value")
        }
    }

}
