/*
 * Copyright (c) 2012, the authors.
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

package nextflow.script

import groovy.util.logging.Slf4j

/**
 * Custom binding used to hold the CLI specified parameters.
 * <p>
 * The difference respect the default implementation is that
 * once the value is defined it cannot be overridden, so this make
 * the parameters definition works like constant values.
 * <p>
 * The main reason for that is to be able to provide optional default value
 * for script parameters in the pipeline script.
 *
 * Read more about 'binding variables'
 * http://groovy.codehaus.org/Scoping+and+the+Semantics+of+%22def%22
 *
 */
@Slf4j
class CliBinding extends Binding {

    def parameters = []

    def void setParam( String name, String value ) {

        // mark this name as a parameter
        if( !parameters.contains(name) ) {
            parameters.add(name)
        }

        super.setVariable(name,parseValue(value))
    }

    def void setVariable(String name, Object value) {

        // variable name marked as parameter cannot be overridden
        if( name in parameters ) {
            return
        }

        super.setVariable(name,value)
    }



    static private parseValue( String str ) {

        if ( str == null ) return null

        if ( str?.toLowerCase() == 'true') return Boolean.TRUE
        if ( str?.toLowerCase() == 'false' ) return Boolean.FALSE

        if ( str.isInteger() ) return str.toInteger()
        if ( str.isLong() ) return str.toLong()
        if ( str.isDouble() ) return str.toDouble()

        return str

    }
}
