/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.k8s.client

import nextflow.util.Duration

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

    Integer maxErrorRetry = 4

    /**
     * Timeout when reading from Input stream when a connection is established to a resource.
     * If the timeout expires before there is data available for read, a {@link java.net.SocketTimeoutException} is raised
     */
    Duration httpReadTimeout

    /**
     * Timeout when opening a communications link to the resource referenced by K8sClient request connection
     * If the timeout expires before there is data available for read, a {@link java.net.SocketTimeoutException} is raised
     */
    Duration httpConnectTimeout

    /**
     * When true signal that the configuration was retrieved from within a K8s cluster
     */
    boolean isFromCluster

    String getNamespace() { namespace ?: 'default' }

    ClientConfig() {

    }

    String toString() {
        "${this.class.getSimpleName()}[ server=$server, namespace=$namespace, serviceAccount=$serviceAccount, token=${cut(token)}, sslCert=${cut(sslCert)}, clientCert=${cut(clientCert)}, clientKey=${cut(clientKey)}, verifySsl=$verifySsl, fromFile=$isFromCluster, httpReadTimeout=$httpReadTimeout, httpConnectTimeout=$httpConnectTimeout, maxErrorRetry=$maxErrorRetry ]"
    }

    private String cut(String str) {
        if( !str ) return '-'
        return str.size()<10 ? str : str[0..10].toString() + '..'
    }

    private String cut(byte[] bytes) {
        if( !bytes ) return '-'
        cut(bytes.encodeBase64().toString())
    }

    static ClientConfig discover(String context, String namespace, String serviceAccount) {
        new ConfigDiscovery().discover(context, namespace, serviceAccount)
    }

    static ClientConfig fromNextflowConfig(Map opts, String namespace, String serviceAccount) {
        final result = new ClientConfig()

        if( opts.server )
            result.server = opts.server

        if( opts.token )
            result.token = opts.token
        else if( opts.tokenFile )
            result.token = Paths.get(opts.tokenFile.toString()).getText('UTF-8')

        result.namespace = namespace ?: opts.namespace ?: 'default'

        result.serviceAccount = serviceAccount ?: 'default'

        if( opts.verifySsl )
            result.verifySsl = opts.verifySsl as boolean

        if( opts.sslCert )
            result.sslCert = opts.sslCert.toString().decodeBase64()
        else if( opts.sslCertFile )
            result.sslCert = Paths.get(opts.sslCertFile.toString()).bytes

        if( opts.clientCert )
            result.clientCert = opts.clientCert.toString().decodeBase64()
        else if( opts.clientCertFile )
            result.clientCert = Paths.get(opts.clientCertFile.toString()).bytes

        if( opts.clientKey )
            result.clientKey = opts.clientKey.toString().decodeBase64()
        else if( opts.clientKeyFile )
            result.clientKey = Paths.get(opts.clientKeyFile.toString()).bytes

        if( opts.maxErrorRetry )
            result.maxErrorRetry = opts.maxErrorRetry as Integer

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
