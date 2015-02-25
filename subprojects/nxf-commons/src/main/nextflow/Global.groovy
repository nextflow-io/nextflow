/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

package nextflow

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Hold global variables
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class Global {

    /**
     * The pipeline session instance
     */
    static private session

    /**
     * The main configuration object
     */
    static Map config

    /**
     * @return The object instance representing the current session
     */
    static <T> T getSession() {
        (T)session
    }

    /**
     * Set the application session object
     *
     * @param value An object instance representing the current session
     */
    static <T> void setSession( value ) {
        session = value
    }

    /**
     * Run the specified closure on application termination
     *
     * @param callback A closure to be executed on application shutdown
     */
    static void onShutdown(Closure callback) {
        hooks.add(callback)
    }

    static final List<Closure> hooks = []

    static private synchronized cleanUp() {
        for( Closure c : hooks ) {
            try {
                c.call()
            }
            catch( Exception e ) {
                log.debug("Error during on cleanup", e )
            }
        }
    }

    /*
     * Global shutdown hook
     */
    static {
        Runtime.getRuntime().addShutdownHook {
            cleanUp()
        }
    }


}
