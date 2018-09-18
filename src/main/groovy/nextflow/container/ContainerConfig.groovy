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

package nextflow.container

import groovy.transform.CompileStatic

/**
 * Models container engine configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ContainerConfig extends LinkedHashMap {

    /* required by Kryo deserialization -- do not remove */
    private ContainerConfig() { }

    ContainerConfig(Map config) {
        super(config)
    }

    boolean isEnabled() {
        get('enabled')?.toString() == 'true'
    }

    boolean isLegacy() {
        get(legacy)?.toString() == 'true'
    }

    String getEngine() {
        get('engine')
    }

    List<String> getEnvWhitelist() {
        def result = get('envWhitelist')
        if( !result )
            return Collections.emptyList()

        if( result instanceof CharSequence )
            return result.tokenize(',').collect { it.trim() }

        if( result instanceof List )
            return result

        throw new IllegalArgumentException("Not a valid `envWhitelist` argument")
    }
}
