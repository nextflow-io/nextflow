/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.util

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import org.apache.commons.lang.StringUtils
/**
 * Helper class retrieve cluster configuration properties
 */
@Slf4j
class ClusterConfig  {

    private String scope

    private Map config

    private Map environment = System.getenv()


    ClusterConfig( Map config, String scope = null, Map env = null ) {
        this.config = config
        this.scope = scope
        this.environment = env
    }

    String getScope() { scope }

    Map getConfig() { config }

    Map getEnvironment() { environment }

    def <T> T getAttribute( String name ) {
        getAttribute(name, null)
    }

    def <T> T getAttribute( String name, defValue ) {
        def result
        if( scope && config instanceof Map && config[scope] instanceof Map ) {
            result = getAttributeValue0( "${scope}.$name", config )
        }

        if( result != null )
            return (T)result

        result = getAttributeValue0( name, config )
        if( result != null )
            return (T)result

        // -- try to fallback on env
        if( environment ) {
            def key = "NXF_CLUSTER_${name.toUpperCase().replaceAll(/\./,'_')}".toString()
            result = environment.get(key)
        }

        return (result != null ? result : defValue) as T
    }

    private Object getAttributeValue0( String key, map ) {

        if( !(map instanceof Map) )
            return null

        def p = key.indexOf('.')
        if( p == -1 )
            return ConfigHelper.parseValue(map.get(key))

        String parent = key.substring(0,p)
        return getAttributeValue0( key.substring(p+1), map.get(parent) )
    }

    Set<String> getAttributeNames() {
        def result = new HashSet(config.keySet())
        def allScopes = result.findAll { String it -> it.startsWith('$') }

        if( allScopes ) {
            result.removeAll(allScopes)
            if( scope ) allScopes = ['$'+scope]
            allScopes.each { key ->
                def map = (Map)config.get(key)
                if( map ) result.addAll( map.keySet() )
            }
        }

        return result
    }

    Set<String> getAttributeNames( String key ) {
        getAttributeNames0(key, config)
    }

    private Set<String> getAttributeNames0( String key, map ) {
        assert key, 'Missing `key` parameter'

        def p = key.indexOf('.')
        if( p != -1 ) {
            String parent = key.substring(0,p)
            return getAttributeNames0( key.substring(p+1), map.get(parent) )
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

    String getClusterJoin() {

        def result = getAttribute('join')

        def seed = environment?.NXF_CLUSTER_SEED as String
        if( !result && seed ) {
            log.debug "Cluster seed number: ${seed}"
            result = seedToMulticastAddress(seed)
        }

        return result
    }

    String getCloudDriverName() {
        if( !isCloudCluster() ) return null
        try {
            // this string is supposed to be in the form: `cloud:<provider>:<cluster-name>`
            return getClusterJoin().tokenize(':')[1]
        }
        catch (Throwable e) {
            log.debug "Oops.. Cannot fetch cloud driver name", e
            return null
        }
    }

    boolean isCloudCluster() {
        getClusterJoin()?.startsWith('cloud:')
    }

    @PackageScope
    static String seedToMulticastAddress( String seed ) {
        try {
            int value = Integer.parseInt(seed)
            if( value<0 || value>= 16777216 ) {
                log.debug "WARN: cluster seed number out of range [0-16777215]: $value"
                return null
            }

            int a = value & 255
            value = value >>> 8
            int b = value & 255
            value = value >>> 8
            int c = value & 255

            return "multicast:228.${c}.${b}.${a}"
        }
        catch( NumberFormatException e ) {
            log.debug "WARN: Not a valid cluster seed number: $seed"
            return null
        }
    }

}