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

package nextflow.executor.res

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import nextflow.util.MemoryUnit

/**
 * Models disk resource request with support for cloud-specific options.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@CompileStatic
@EqualsAndHashCode
class DiskResource {

    final MemoryUnit request
    final String type
    final Integer iops
    final Integer throughput
    final Boolean encrypted
    final String filesystem
    final String mountPath

    DiskResource( value ) {
        this(request: value)
    }

    DiskResource( Map opts ) {
        this.request = toMemoryUnit(opts.request)

        if( opts.type )
            this.type = opts.type as String
        if( opts.iops )
            this.iops = opts.iops as Integer
        if( opts.throughput )
            this.throughput = opts.throughput as Integer
        if( opts.encrypted != null )
            this.encrypted = opts.encrypted as Boolean
        if( opts.filesystem )
            this.filesystem = opts.filesystem as String
        if( opts.mountPath )
            this.mountPath = opts.mountPath as String
    }

    DiskResource withRequest(MemoryUnit value) {
        return new DiskResource(
            request: value,
            type: this.type,
            iops: this.iops,
            throughput: this.throughput,
            encrypted: this.encrypted,
            filesystem: this.filesystem,
            mountPath: this.mountPath
        )
    }

    private static MemoryUnit toMemoryUnit( value ) {
        if( value instanceof MemoryUnit )
            return (MemoryUnit)value

        try {
            return new MemoryUnit(value.toString().trim())
        }
        catch( Exception e ) {
            throw new IllegalArgumentException("Not a valid disk value: $value")
        }
    }

}
