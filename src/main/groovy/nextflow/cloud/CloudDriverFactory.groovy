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

package nextflow.cloud

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.util.ServiceDiscover
import nextflow.util.ServiceName

/**
 * Factory class to get and instance of a concrete {@link CloudDriver} object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CloudDriverFactory {

    /**
     * @return The map of cloud drivers available on the current class path
     */
    @Memoized
    static Map<String,Class<CloudDriver>> loadDrivers() {

        def result = new HashMap<String,Class<CloudDriver>>()
        ServiceDiscover.load(CloudDriver).iterator().each { Class<CloudDriver> driver ->
            final name = nameForClass(driver)
            log.debug "Discovered cloud driver: `$name` [${driver.getName()}]"
            result.put(name, driver)
        }

        return result
    }

    private static String nameForClass(Class driverClass) {
        def name = driverClass.getAnnotation(ServiceName)
        return name ? name.value() : driverClass.getSimpleName().toLowerCase().replace('driver','')
    }

    /**
     * @return The list of name of available drivers
     */
    static Set<String> getDriverNames() {
        return loadDrivers().keySet()
    }

    /**
     * @return The first and only cloud driver name or {@code null} otherwise
     */
    static String getDefaultDriverName() {
        def names = getDriverNames()
        names.size() == 1 ? names.iterator().next() : null
    }

    /**
     * Get and instance of a concrete {@link CloudDriver} by the given driver name
     *
     * @param name The name of the cloud driver e.g. {@code aws}
     * @return The cloud driver instance for the given name
     * @throws IllegalArgumentException if the driven with the specified name does not exist
     */
    static CloudDriver getDriver(String name, Map config=null) {
        def result = loadDrivers().get(name)
        if( !result ) throw new IllegalArgumentException("Unknown cloud driver name: `$name`")
        return config ? result.newInstance(config) : result.newInstance()
    }

}
