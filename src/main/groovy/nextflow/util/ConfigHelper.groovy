/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.util
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 * Helper method to handle configuration object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ConfigHelper {


    def static getConfigProperty( def config, String execName, String propName ) {
        def result = null

        // make sure that the *executor* is a map object
        // it could also be a plain string (when it specifies just the its name)
        if( execName && config instanceof Map && config['$'+execName] instanceof Map ) {
            result = config['$'+execName][propName]
        }

        if( result==null && config instanceof Map && config[propName] ) {
            result = config[propName]
        }

        return result
    }

    /**
     * Given a string value converts to its native object representation.
     *
     * @param str A string that may represent boolean, integer, {@link Duration} values
     * @return A object representing the argument value of the string itself if it was not a boolean/number/duration value
     */
    static parseValue( String str ) {

        if ( str == null ) return null

        if ( str.toLowerCase() == 'true') return Boolean.TRUE
        if ( str.toLowerCase() == 'false' ) return Boolean.FALSE

        if ( str.isInteger() ) return str.toInteger()
        if ( str.isLong() ) return str.toLong()
        if ( str.isDouble() ) return str.toDouble()
        // try as duration as well
        try { return new Duration(str) }
        catch( IllegalArgumentException e ) { }

        return str

    }

    static parseValue( obj ) {
        if( obj instanceof String )
            return parseValue((String)obj)

        if( obj instanceof GString )
            return parseValue(obj.toString())

        return obj
    }

    /**
     * Given a list of paths looks for the files ending withe the extension '.jar' and return
     * a list containing the original directories, plus the JARs paths
     *
     * @param dirs
     * @return
     */
    static List<Path> resolveClassPaths( List<Path> dirs ) {

        List<Path> result = []
        if( !dirs )
            return result

        for( Path path : dirs ) {
            if( path.isFile() && path.name.endsWith('.jar') ) {
                result << path
            }
            else if( path.isDirectory() ) {
                result << path
                path.eachFileMatch( ~/.+\.jar$/ ) { if(it.isFile()) result << it }
            }
        }

        return result
    }



}

