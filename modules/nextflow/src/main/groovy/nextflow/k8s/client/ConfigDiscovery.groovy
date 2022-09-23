/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.k8s.client

import javax.net.ssl.KeyManager
import javax.net.ssl.KeyManagerFactory
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.security.KeyStore

import groovy.util.logging.Slf4j
import org.yaml.snakeyaml.Yaml
/**
 * Discover Kubernetes configuration from system environment
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ConfigDiscovery {

    static final ConfigDiscovery instance = new ConfigDiscovery()

    static ConfigDiscovery getInstance() { instance }

    /**
     * Discover Kubernetes service from current environment settings
     */
    ClientConfig discover(String context=null) {

        // Note: System.getProperty('user.home') may not report the correct home path when
        // running in a container. Use env HOME instead.
        def env = System.getenv()
        def kubeConfig = env.get('KUBECONFIG') ?: "${env.get('HOME')}/.kube/config"
        def configFile = Paths.get(kubeConfig)

        // discover client config from KUBECONFIG (i.e. running in local machine)
        if( configFile.exists() ) {
            return fromConfig(configFile, context)
        }
        else {
            log.debug "K8s config file does not exist: $configFile"
        }

        // discover client config from K8s env vars (i.e. running in a pod)
        if( env.get('KUBERNETES_SERVICE_HOST') ) {
            return fromCluster(env)
        }
        else {
            log.debug "K8s env variable KUBERNETES_SERVICE_HOST is not defined"
        }

        throw new IllegalStateException("Unable to lookup Kubernetes cluster configuration")
    }

    protected ClientConfig fromCluster(Map<String,String> env) {

        // See https://kubernetes.io/docs/tasks/run-application/access-api-from-pod/

        final host = env.get('KUBERNETES_SERVICE_HOST')
        final port = env.get('KUBERNETES_SERVICE_PORT')
        final server = host + ( port ? ":$port" : '' )

        final cert = path('/var/run/secrets/kubernetes.io/serviceaccount/ca.crt').bytes
        final token = path('/var/run/secrets/kubernetes.io/serviceaccount/token').text
        final namespace = path('/var/run/secrets/kubernetes.io/serviceaccount/namespace').text

        new ClientConfig( server: server, token: token, namespace: namespace, sslCert: cert, isFromCluster: true )
    }

    /**
     * Helper method for testing purposes.
     */
    protected Path path(String path) {
        Paths.get(path)
    }

    protected ClientConfig fromConfig(Path path, String contextName=null) {
        def yaml = (Map)new Yaml().load(Files.newInputStream(path))

        contextName = contextName ?: yaml."current-context" as String

        final allContexts = yaml.contexts as List
        final allClusters = yaml.clusters as List
        final allUsers = yaml.users as List
        final context = allContexts.find { Map it -> it.name == contextName } ?.context
        if( !context )
            throw new IllegalArgumentException("Unknown Kubernetes context: $contextName -- check config file: $path")
        final userName = context?.user
        final clusterName = context?.cluster
        final user = allUsers.find{ Map it -> it.name == userName } ?.user ?: [:]
        final cluster = allClusters.find{ Map it -> it.name == clusterName } ?.cluster ?: [:]

        def config = ClientConfig.fromUserAndCluster(user, cluster, path)

        if( context?.namespace ) {
            config.namespace = context?.namespace
        }

        return config
    }

    protected KeyStore createKeyStore0(byte[] clientCert, byte[] clientKey, char[] passphrase, String alg) {
        def cert = new ByteArrayInputStream(clientCert)
        def key = new ByteArrayInputStream(clientKey)
        return SSLUtils.createKeyStore(cert, key, alg, passphrase, null, null)
    }

    protected KeyStore createKeyStore(byte[] clientCert, byte[] clientKey, char[] passphrase) {
        try {
            // try first RSA algorithm
            return createKeyStore0(clientCert, clientKey, passphrase, "RSA")
        }
        catch (Exception e1) {
            // fallback to EC algorithm
            try {
                return createKeyStore0(clientCert, clientKey, passphrase, "EC")
            }
            catch (Exception e2) {
                // if still fails, throws the first exception
                throw e1
            }
        }
    }

    KeyManager[] createKeyManagers(byte[] clientCert, byte[] clientKey) {
        final passphrase = "".toCharArray()
        final keyStore = createKeyStore(clientCert, clientKey, passphrase)
        final kmf = KeyManagerFactory.getInstance(KeyManagerFactory.getDefaultAlgorithm());
        kmf.init(keyStore, passphrase);
        return kmf.getKeyManagers();
    }

    String discoverAuthToken(String namespace, String serviceAccount) {
        def cmd = "kubectl -n ${namespace} get secret `kubectl -n ${namespace} get serviceaccount ${serviceAccount} -o jsonpath='{.secrets[0].name}'` -o jsonpath='{.data.token}'"
        def proc = new ProcessBuilder('bash','-o','pipefail','-c', cmd).redirectErrorStream(true).start()
        def status = proc.waitFor()
        if( status==0 ) {
            return new String(proc.text.trim().decodeBase64())
        }
        else {
            final msg = proc.text
            final cause = msg ? "-- cause:\n${msg.indent('  ')}" : ''
            log.debug "[K8s] unable to fetch auth token ${cause}"
            return null
        }
    }

}
