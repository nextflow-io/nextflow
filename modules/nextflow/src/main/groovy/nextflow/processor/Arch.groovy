/*
 * Copyright 2022, Pawsey Supercomputing Research Centre
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

package nextflow.processor

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Implements the {@code arch} process directive, to hold information on the 
 * CPU (micro)architecture required by the process.
 * 
 * @author Marco De La Pierre <marco.delapierre@gmail.com>
 */
@ToString
@EqualsAndHashCode
@CompileStatic
class Arch {

/*
 * example of notation in process: arch 'linux/x86_64', target: 'haswell'
 * example of notation in config:  arch = [name: 'linux/x86_64', target: 'haswell']
 * 
 * where fullArch = 'linux/x86_64'
 *       platform = 'linux'
 *       arch = 'x86_64'
 *       target = 'haswell'
 * 
 *       spackArch = target ?: arch
 */
    final String fullArch
    final String platform
    final String arch
    final String target
    final String spackArch

    @CompileStatic
    protected String getPlatform( String value ) {
        return value.minus(~'/.*')
    }

    @CompileStatic
    protected String getArch( String value ) {
        return value.minus(~'.*/')
    }

    @CompileStatic
    protected String getSpackArch( Map res ) {
        if( res.target != null )
            return res.target as String
        else if( res.name != null )
            return getArch(res.name as String)
        else
            return null
    }

    Arch( String value ) {
        this(name: value)
    }

    Arch( Map res ) {
        if( res.name != null ) {
            this.fullArch = res.name as String
            this.platform = getPlatform(res.name as String)
            this.arch = getArch(res.name as String)
        }
        if( res.target != null )
            this.target = res.target as String
        if( res.name!=null || res.target!=null )
            this.spackArch = getSpackArch(res)
    }

}
