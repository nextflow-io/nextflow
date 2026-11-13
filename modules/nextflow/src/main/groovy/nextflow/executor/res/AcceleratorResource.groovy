/*
 * Copyright 2013-2026, Seqera Labs
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
import java.util.regex.Pattern

/**
 * Models accelerator resource request
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@CompileStatic
@EqualsAndHashCode
class AcceleratorResource {

    private static final Pattern CUDA_VERSION = ~/^\d+\.\d+$/

    final Integer request
    final Integer limit
    final String type
    final String runtime

    /**
     * The CUDA release (major.minor, e.g. {@code 12.2}) the task container image was built with
     */
    final String cuda

    AcceleratorResource( Number value ) {
        this(limit: value)
    }

    AcceleratorResource( Map res ) {
        if( res.limit != null && res.request != null ) {
            this.limit = res.limit as int
            this.request = res.request as int
        }
        else if( res.limit != null ) {
            this.limit = res.limit as int
            this.request = res.limit as int
        }
        else if( res.request != null ) {
            this.request = res.request as int
        }
        else {
            throw new IllegalArgumentException("Invalid `accelerator` directive value -- `request` or `limit` is required: $res")
        }

        if( res.type )
            this.type = res.type as String

        if( res.runtime )
            this.runtime = res.runtime as String

        if( res.cuda != null )
            this.cuda = parseCuda(res.cuda)
    }

    static protected String parseCuda(Object value) {
        // a number is rejected because it would be lossy, e.g. `12.10` would become `12.1`
        if( value !instanceof CharSequence )
            throw new IllegalArgumentException("Invalid `accelerator` directive `cuda` option -- it must be a string in the form 'major.minor' e.g. cuda: '12.2' (offending value: $value)")
        final result = value.toString().trim()
        if( !CUDA_VERSION.matcher(result).matches() )
            throw new IllegalArgumentException("Invalid `accelerator` directive `cuda` option -- it must be a CUDA release in the form 'major.minor' e.g. cuda: '12.2' (offending value: '$value')")
        return result
    }

}
