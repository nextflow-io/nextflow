/*
 * Copyright 2023, Pawsey Supercomputing Research Centre
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
class Architecture {

/*
 * example of notation in process: arch 'linux/x86_64', target: 'haswell'
 * example of notation in config:  arch = [name: 'linux/x86_64', target: 'haswell']
 * 
 * where dockerArch = 'linux/x86_64'
 *       spackArch = target ?: arch  // plus some validation for Spack syntax
 * 
 *       platform = 'linux'
 *       arch = 'x86_64'
 *       target = 'haswell'
 * 
 * [alternate example: 'arch linux/arm/v8', where platform = 'linux' and arch = 'arm/v8']
 */
    // used in Nextflow
    final String dockerArch
    final String spackArch

    // defined, but currently not used
    final String platform
    final String arch
    final String target

    static protected String getPlatform( String value ) {
        // return value.minus(~'/.*') // keeping for reference
        final chunks = value.tokenize('/')
        if( chunks.size() > 1 )
            return chunks[0]
        else
            return null
    }

    static protected String getArch( String value ) {
        // return value.minus(~'.*/') // keeping for reference
        def chunks = value.tokenize('/')
        if( chunks.size() == 3 )
            return chunks[1] + '/' + chunks[2]
        else if( chunks.size() == 2 )
            return chunks[1]
        else
            return chunks[0]
    }

    static private String validateArchToDockerArch( Map res ) {
        def value = getArch(res.name as String)
        def name = res.name as String
        if( value == 'x86_64' || value == 'amd64' )
            return 'linux/amd64'
        if( value == 'aarch64' || value == 'arm64' || value == 'arm64/v8' )
            return 'linux/arm64'
        if( value == 'arm64/v7' )
            return 'linux/arm64/v7'
        throw new IllegalArgumentException("Not a valid `arch` value: ${name}")
    }

    static private String validateArchToSpackArch( String value, String inputArch ) {
        if( value == 'x86_64' || value == 'amd64' )
            return 'x86_64'
        if( value == 'aarch64' || value == 'arm64' || value == 'arm64/v8' )
            return 'aarch64'
        if( value == 'arm64/v7' )
            return null
        throw new IllegalArgumentException("Not a valid `arch` value: ${inputArch}")
    }

    static protected String getSpackArch( Map res ) {
        if( res.target != null )
            return res.target as String
        else if( res.name != null )
            return validateArchToSpackArch(getArch(res.name as String), res.name as String)
        else
            return null
    }

    Architecture( String value ) {
        this(Map.of('name', value))
    }

    Architecture( Map res ) {
        if( !res.name )
            throw new IllegalArgumentException("Missing architecture `name` attribute")

        this.dockerArch = validateArchToDockerArch(res)
        this.platform = getPlatform(res.name as String)
        this.arch = getArch(res.name as String)

        if( res.target != null )
            this.target = res.target as String
        if( res.name!=null || res.target!=null )
            this.spackArch = getSpackArch(res)
    }

}
