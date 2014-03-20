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

package nextflow.util

import groovy.util.logging.Slf4j

/**
 * Helper method to handle configuration object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ConfigHelper {


    def static getConfigProperty( Map config, String execName, String propName, defValue ) {
        def result = null

        // make sure that the *executor* is a map object
        // it could also be a plain string (when it specifies just the its name)
        if( execName && config instanceof Map && config['$'+execName] instanceof Map ) {
            result = config['$'+execName][propName]
        }

        if( result==null && config instanceof Map && config[propName] ) {
            result = config[propName]
        }


        if( result==null ) {
            result = defValue
            log.trace "Undefined executor property: '$propName' -- fallback default value: $result"
        }

        return result
    }

}
