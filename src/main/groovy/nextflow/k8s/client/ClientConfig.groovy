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

package nextflow.k8s.client
import javax.net.ssl.KeyManager
import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
/**
 * Models the kubernetes cluster client configuration settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@EqualsAndHashCode
@CompileStatic
class ClientConfig {

    boolean verifySsl

    String server

    String namespace

    /**
     * k8s service account name
     * https://kubernetes.io/docs/tasks/configure-pod-container/configure-service-account/
     */
    String serviceAccount

    String token

    byte[] sslCert

    byte[] clientCert

    byte[] clientKey

    KeyManager[] keyManagers

    /**
     * When true signal that the configuration was retrieved from within a K8s cluster
     */
    boolean isFromCluster

    String getNamespace() { namespace ?: 'default' }

    ClientConfig() {

    }

    String toString() {
        "${this.class.getSimpleName()}[ server=$server, namespace=$namespace, token=${cut(token)}, sslCert=${cut(sslCert)}, clientCert=${cut(clientCert)}, clientKey=${cut(clientKey)}, verifySsl=$verifySsl, fromFile=$isFromCluster ]"
    }

    private String cut(String str) {
        if( !str ) return '-'
        return str.size()<10 ? str : str[0..10].toString() + '..'
    }

    private String cut(byte[] bytes) {
        if( !bytes ) return '-'
        cut(bytes.encodeBase64().toString())
    }

    static ClientConfig discover(String context=null) {
        new ConfigDiscovery(context: context).discover()
    }

    static ClientConfig fromMap(Map map) {
        def result = new ClientConfig()
        if( map.server )
            result.server = map.server

        if( map.token )
            result.token = map.token
        else if( map.tokenFile )
            result.token = Paths.get(map.tokenFile.toString()).getText('UTF-8')

        if( map.namespace )
            result.namespace = map.namespace

        if( map.verifySsl )
            result.verifySsl = map.verifySsl as boolean

        if( map.sslCert )
            result.sslCert = map.sslCert.toString().decodeBase64()
        else if( map.sslCertFile )
            result.sslCert = Paths.get(map.sslCertFile.toString()).bytes

        if( map.clientCert )
            result.clientCert = map.clientCert.toString().decodeBase64()
        else if( map.clientCertFile )
            result.clientCert = Paths.get(map.clientCertFile.toString()).bytes

        if( map.clientKey )
            result.clientKey = map.clientKey.toString().decodeBase64()
        else if( map.clientKeyFile )
            result.clientKey = Paths.get(map.clientKeyFile.toString()).bytes

        return result
    }

    static ClientConfig fromUserAndCluster(Map user, Map cluster, Path location) {
        final base = location.isDirectory() ? location : location.parent
        final result = new ClientConfig()
        if( user.token )
            result.token = user.token

        else if( user.tokenFile ) {
            result.token = Paths.get(user.tokenFile.toString()).getText('UTF-8')
        }

        if( user."client-certificate" )
            result.clientCert = base.resolve(user."client-certificate".toString()).bytes

        else if( user."client-certificate-data" )
            result.clientCert = user."client-certificate-data".toString().decodeBase64()

        if( user."client-key" )
            result.clientKey = base.resolve(user."client-key".toString()).bytes

        else if( user."client-key-data" )
            result.clientKey = user."client-key-data".toString().decodeBase64()

        // -- cluster settings

        if( cluster.server )
            result.server = cluster.server

        if( cluster."certificate-authority-data" )
            result.sslCert = cluster."certificate-authority-data".toString().decodeBase64()

        else if( cluster."certificate-authority" )
            result.sslCert = base.resolve(cluster."certificate-authority".toString()).bytes

        result.verifySsl = cluster."insecure-skip-tls-verify" != true

        return result
    }

}
