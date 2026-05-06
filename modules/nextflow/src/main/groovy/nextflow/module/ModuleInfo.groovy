/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.module

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

import java.nio.file.Files
import java.nio.file.Path

/**
 * Utility class for managing .module-info
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class ModuleInfo {

    public static final String MODULE_INFO_FILE = ".module-info"


    /**
     * Save a property to the .module-info file in the module directory
     *
     * @param moduleDir The module directory path
     * @param property The property to save
     * @param value The property value to save
     */
    static void save(Path moduleDir, String property, String value) {
        def moduleInfoFile = moduleDir.resolve(MODULE_INFO_FILE)
        def props = new Properties()
        // If file exists loads to update current just checksum property
        if( Files.exists( moduleInfoFile))
            moduleInfoFile.withInputStream { is -> props.load(is) }
        props.setProperty(property, value)
        moduleInfoFile.withOutputStream { os -> props.store(os, null) }
    }

    /**
     * Save a property to the .module-info file in the module directory
     *
     * @param moduleDir The module directory path
     * @param properties Map with properties to save
     */
    static void save(Path moduleDir, Map<String,String> properties) {
        if( properties ) {
            def moduleInfoFile = moduleDir.resolve(MODULE_INFO_FILE)
            def props = new Properties()
            // If file exists loads to update current just checksum property
            if( Files.exists( moduleInfoFile ) )
                moduleInfoFile.withInputStream { is -> props.load(is) }

            for( final property : properties.entrySet() ) {
                props.setProperty(property.key, property.value)
            }
            moduleInfoFile.withOutputStream { os -> props.store(os, null) }
        }
    }

    /**
     * Return the value of property from the .module-info file in the module directory
     *
     * @param moduleDir The module directory path
     * @param moduleDir The module directory path
     * @return The checksum, or null if file doesn't exist
     */
    static String load(Path moduleDir, String property) {
        def moduleInfoFile = moduleDir.resolve(MODULE_INFO_FILE)
        if( !Files.exists(moduleInfoFile) ) {
            log.debug("Module file $moduleInfoFile not found")
            return null
        }
        def props = new Properties()
        moduleInfoFile.withInputStream { is -> props.load(is) }
        return props.getProperty(property)
    }

    /**
     * Load all properties from the .module-info file in the module directory
     *
     * @param moduleDir The module directory path
     * @return Map of properties, or null if file doesn't exist
     */
    static Map<String, String> load(Path moduleDir) {
        def moduleInfoFile = moduleDir.resolve(MODULE_INFO_FILE)
        if( !Files.exists(moduleInfoFile) ) {
            log.debug("Module file $moduleInfoFile not found")
            return null
        }
        def props = new Properties()
        moduleInfoFile.withInputStream { is -> props.load(is) }
        return props as Map<String, String>
    }
}
