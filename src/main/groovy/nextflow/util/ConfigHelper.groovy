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

package nextflow.util

import groovy.util.logging.Slf4j
import org.apache.commons.lang.StringUtils

/**
 * Helper method to handle configuration object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
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

}

/**
 * Helper class retrieve daemon configuration properties
 */
class DaemonConfig {

    final String scope

    final Map config

    final Map fallbackEnvironment


    DaemonConfig( String scope, Map config, Map env = null ) {
        assert scope
        assert config != null

        this.scope = scope
        this.config = config
        this.fallbackEnvironment = env
    }

    def <T> T getAttribute( String name, defValue = null ) {
        def result = ConfigHelper.getConfigProperty(config, scope, name)
        if( result != null )
            return result

        // -- try to fallback on env
        if( fallbackEnvironment ) {
            def key = "NXF_DAEMON_${name.toUpperCase().replaceAll(/\./,'_')}".toString()
            result = fallbackEnvironment.get(key)
        }

        return result != null ? result : defValue
    }



    /**
     * Get the network interface IP addresses. When the it specified by its name
     * it will resolve to the actual IP address
     *
     * @return The list of address or an empty list if attribute is not specified
     */
    List<String> getNetworkInterfaceAddresses() {
        def result = []
        def value = getAttribute('interface') as String
        def interfaceNames = value ? StringUtils.split(value, ", \n").collect { it.trim() } : null
        interfaceNames?.each {
            if( it.contains('.') )
                result << it // it is supposed to be an interface IP address, add it to the list

            else
                result.addAll( findInterfaceAddressesByName(it) )
        }
        return result
    }


    protected List<String> findInterfaceAddressesByName( String name ) {
        assert name
        def result = []
        NetworkInterface interfaces = NetworkInterface.getNetworkInterfaces().find { NetworkInterface net -> net.name == name || net.displayName == name }
        interfaces?.getInetAddresses()?.each { InetAddress addr ->
            if( addr instanceof Inet4Address )
                result << addr.getHostAddress()
        }
        return result
    }


}