/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
        ServiceDiscover.load(CloudDriver).each { Class<CloudDriver> driver ->
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
