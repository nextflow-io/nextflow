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

import javax.net.ssl.KeyManager
import javax.net.ssl.KeyManagerFactory
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.security.KeyStore
import static nextflow.util.StringUtils.formatHostName

import groovy.util.logging.Slf4j
import org.yaml.snakeyaml.Yaml
/**
 * Discover Kubernetes configuration from system environment
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ConfigDiscovery {

    private Map<String,String> env = System.getenv()

    ConfigDiscovery() { }

    /**
     * Discover Kubernetes client configuration from current environment using
     * either the .kube/config file or the pod service service account virtual
     * file system when running in a pod.
     *
     * @param contextName The K8s config context name.
     * @param namespace The K8s cluster namespace
     * @param serviceAccount The K8s cluster service account
     * @return The e
     */
    ClientConfig discover(String contextName, String namespace, String serviceAccount) {

        // Note: System.getProperty('user.home') may not report the correct home path when
        // running in a container. Use env HOME instead.
        def home = System.getenv('HOME')
        def kubeConfig = env.get('KUBECONFIG') ? env.get('KUBECONFIG') : "$home/.kube/config"
        def configFile = Paths.get(kubeConfig)

        // determine the Kubernetes client configuration via the `.kube/config` file
        if( configFile.exists() ) {
            return fromKubeConfig(configFile, contextName, namespace, serviceAccount)
        }
        else {
            log.debug "K8s config file does not exist: $configFile"
        }

        // determine the Kubernetes client configuration via the pod environment
        if( env.get('KUBERNETES_SERVICE_HOST') ) {
            return fromCluster(env, namespace, serviceAccount)
        }
        else {
            log.debug "K8s env variable KUBERNETES_SERVICE_HOST is not defined"
        }

        throw new IllegalStateException("Unable to lookup Kubernetes cluster configuration")
    }

    protected ClientConfig fromCluster(Map<String,String> env, String cfgNamespace, String serviceAccount) {

        // See https://kubernetes.io/docs/tasks/access-application-cluster/access-cluster/#accessing-the-api-from-a-pod

        final host = env.get('KUBERNETES_SERVICE_HOST')
        final port = env.get('KUBERNETES_SERVICE_PORT')
        final server = formatHostName(host, port)

        final cert = path('/var/run/secrets/kubernetes.io/serviceaccount/ca.crt').bytes
        final token = path('/var/run/secrets/kubernetes.io/serviceaccount/token').text
        final namespace = path('/var/run/secrets/kubernetes.io/serviceaccount/namespace').text

        return new ClientConfig(
                server: server,
                token: token,
                namespace: cfgNamespace ?: namespace,
                serviceAccount: serviceAccount,
                sslCert: cert,
                isFromCluster: true )
    }

    protected Path path(String path) {
       Paths.get(path)
    }

    protected ClientConfig fromKubeConfig(Path path, String contextName, String namespace, String serviceAccount) {
        def yaml = (Map)new Yaml().load(Files.newInputStream(path))

        contextName ?= yaml."current-context" as String

        final allContext = yaml.contexts as List
        final allClusters = yaml.clusters as List
        final allUsers = yaml.users as List
        final context = allContext.find { Map it -> it.name == contextName } ?.context
        if( !context )
            throw new IllegalArgumentException("Unknown Kubernetes context: $contextName -- check config file: $path")
        final userName = context?.user
        final clusterName = context?.cluster
        final user = allUsers.find{ Map it -> it.name == userName } ?.user ?: [:]
        final cluster = allClusters.find{ Map it -> it.name == clusterName } ?.cluster ?: [:]

        final config = ClientConfig.fromUserAndCluster(user, cluster, path)

        // the namespace provided should have priority over the context current namespace
        config.namespace = namespace ?: context?.namespace ?: 'default'

        config.serviceAccount = serviceAccount ?: 'default'

        if( config.clientCert && config.clientKey ) {
            config.keyManagers = createKeyManagers(config.clientCert, config.clientKey)
        }
        else if( !config.token ) {
            config.token = discoverAuthToken(contextName, config.namespace, config.serviceAccount)
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

    protected KeyManager[] createKeyManagers(byte[] clientCert, byte[] clientKey) {
        final passphrase = "".toCharArray()
        final keyStore = createKeyStore(clientCert, clientKey, passphrase)
        final kmf = KeyManagerFactory.getInstance(KeyManagerFactory.getDefaultAlgorithm());
        kmf.init(keyStore, passphrase);
        return kmf.getKeyManagers();
    }

    String discoverAuthToken(String context, String namespace, String serviceAccount) {
        context ?= 'default'
        namespace ?= 'default'
        serviceAccount ?= 'default'

        final cmd = "kubectl --context $context -n ${namespace} get secret -o=jsonpath='{.items[?(@.metadata.annotations.kubernetes\\.io/service-account\\.name==\"$serviceAccount\")].data.token}'"
        final proc = new ProcessBuilder('bash','-o','pipefail','-c', cmd).start()
        final status = proc.waitFor()
        final text = proc.inputStream?.text
        if( status==0 && text ) {
            try {
                return new String(text.trim().decodeBase64())
            }
            catch( Exception e ) {
                log.warn "Unable to decode K8s cluster auth token '$text' -- cause: ${e.message}"
            }
        }
        else {
            final cause = proc.errorStream?.text ?: text
            final msg = cause ? "\n- cmd  : $cmd\n- exit : $status\n- cause:\n${cause.indent('  ')}" : ''
            log.warn "[K8s] unable to fetch auth token ${msg}"
        }
        return null
    }

}
