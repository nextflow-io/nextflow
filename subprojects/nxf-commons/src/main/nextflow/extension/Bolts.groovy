/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.extension

import groovy.transform.CompileStatic

/**
 * Nextflow extensions
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Bolts {

    /**
     * Converts {@code ConfigObject}s to a plain {@code Map}
     *
     * @param config
     * @return A normalized config object
     */
    static Map toMap( ConfigObject config ) {
        assert config != null
        (Map)normalize0((Map)config)
    }

    static private normalize0( config ) {

        if( config instanceof Map ) {
            Map result = new LinkedHashMap(config.size())
            config.keySet().each { name ->
                def value = (config as Map).get(name)
                result.put(name, normalize0(value))
            }
            return result
        }
        else if( config instanceof Collection ) {
            List result = new ArrayList(config.size())
            for( entry in config ) {
                result << normalize0(entry)
            }
            return result
        }
        else if( config instanceof GString ) {
            return config.toString()
        }
        else {
            return config
        }
    }

}
