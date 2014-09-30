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

import java.nio.file.Path

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
        def result
        if( scope && config instanceof Map && config['$'+scope] instanceof Map ) {
            result = getAttrValue( "\$${scope}.$name", config )
        }

        if( result != null )
            return (T)result

        result = getAttrValue( name, config )
        if( result != null )
            return (T)result

        // -- try to fallback on env
        if( fallbackEnvironment ) {
            def key = "NXF_CLUSTER_${name.toUpperCase().replaceAll(/\./,'_')}".toString()
            result = fallbackEnvironment.get(key)
        }

        return (T) result != null ? result : defValue
    }

    private Object getAttrValue( String key, map ) {

        if( !(map instanceof Map) )
            return null

        def p = key.indexOf('.')
        if( p == -1 )
            return ConfigHelper.parseValue(map.get(key))

        String parent = key.substring(0,p)
        return getAttrValue( key.substring(p+1), map.get(parent) )
    }

    def List<String> getAttributesNames( String key = null ) {
        getAttrNames(key, config)
    }

    private List<String> getAttrNames( String key, map ) {

        if( !key ) {
           return map instanceof Map ? new ArrayList<>(map.keySet()) : []
        }

        def p = key.indexOf('.')
        if( p != -1 ) {
            String parent = key.substring(0,p)
            return getAttrNames( key.substring(p+1), map.get(parent) )
        }
        else {
            def entry = map.get(key)
            entry instanceof Map ? new ArrayList<>(entry.keySet()) : []
        }

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