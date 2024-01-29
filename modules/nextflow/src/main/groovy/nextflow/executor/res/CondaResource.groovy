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
 *
 */

package nextflow.executor.res


import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
/**
 * Model Conda packaging tool resource
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
class CondaResource {

    private static final List<String> ALLOWED_KEYS = List.<String>of('packages','pip')

    final String packages
    final String pip

    private CondaResource(String conda, String pip) {
        if( pip && conda && isFilePath(conda) )
            throw new IllegalArgumentException("Conda environment file and Pip packages conflict each other - offending enviroment: '$conda'")
        this.packages = conda
        this.pip = pip
    }

    String getExtendedPackages() {
        if( !pip )
            return packages
        if( !packages )
            return prefix(pip)
        else
            return packages + ' ' + prefix(pip)
    }

    protected String prefix(String value) {
        if( !value )
            return null
        value.tokenize(' ').collect(it-> "pip:"+it).join(' ')
    }

    static CondaResource ofCondaPackages(CharSequence packages) {
        new CondaResource(packages.toString(), null)
    }

    static CondaResource ofPipPackages(CharSequence packages) {
        new CondaResource(null, packages.toString())
    }

    static CondaResource of(Map<String,?> values) {
        if( !values )
            throw new IllegalArgumentException("Conda directive cannot be empty")
        for( String k in values.keySet()) {
            if( k !in ALLOWED_KEYS )
                throw new IllegalArgumentException("Not a valid Conda directive attribute: '$k'")
        }
        return new CondaResource(values.packages?.toString(), values.pip?.toString())
    }

    static boolean isFilePath(String value) {
        if( !value )
            return false
        if( value.contains('/') )
            return true
        if( value.endsWith('.txt') )
            return true
        if( value.endsWith('.yml') )
            return true
        if( value.endsWith('.yaml') )
            return true
        if( value.startsWith('http://'))
            return true
        if( value.startsWith('https://'))
            return true
        return false
    }
}
