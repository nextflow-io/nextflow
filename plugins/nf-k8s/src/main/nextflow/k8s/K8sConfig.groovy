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

package nextflow.k8s

import nextflow.k8s.client.K8sRetryConfig

import java.util.concurrent.TimeUnit
import javax.annotation.Nullable

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.BuildInfo
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.container.ContainerHelper
import nextflow.script.dsl.Description
import nextflow.exception.AbortOperationException
import nextflow.k8s.client.ClientConfig
import nextflow.k8s.client.K8sClient
import nextflow.k8s.client.K8sResponseException
import nextflow.k8s.model.PodOptions
import nextflow.k8s.model.PodSecurityContext
import nextflow.k8s.model.PodVolumeClaim
import nextflow.k8s.model.ResourceType
import nextflow.util.Duration
/**
 * Model Kubernetes specific settings defined in the nextflow
 * configuration file
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("k8s")
@Description("""
    The `k8s` scope controls the deployment and execution of workflow applications in a Kubernetes cluster.
""")
@Slf4j
@CompileStatic
class K8sConfig implements ConfigScope {

    static final private Map<String,?> DEFAULT_FUSE_PLUGIN = Map.of('nextflow.io/fuse', 1)

    private Cache<String, ClientConfig> clientCache

    @ConfigOption
    @Description("""
        Automatically mount host paths into the task pods (default: `false`). Only intended for development purposes when using a single node.
    """)
    final boolean autoMountHostPaths

    @ConfigOption
    @Description("""
        Whether to use Kubernetes `Pod` or `Job` resource type to carry out Nextflow tasks (default: `Pod`).
    """)
    final String computeResourceType

    @ConfigOption
    @Description("""
        When `true`, successful pods are automatically deleted (default: `true`).
    """)
    final private Boolean cleanup

    @ConfigOption
    @Description("""
        Map of options for the K8s client.

        If this option is specified, it will be used instead of `.kube/config`.
    """)
    final Map client

    @ConfigOption
    @Description("""
        The interval after which the Kubernetes client configuration is refreshed (default: `50m`).
    """)
    final Duration clientRefreshInterval

    @ConfigOption
    @Description("""
        The Kubernetes [configuration context](https://kubernetes.io/docs/tasks/access-application-cluster/configure-access-multiple-clusters/) to use.
    """)
    final String context

    @ConfigOption
    @Description("""
        When `true`, set both the pod CPUs `request` and `limit` to the value specified by the `cpus` directive, otherwise set only the `request` (default: `false`).
    """)
    final boolean cpuLimits

    final K8sDebug debug

    @ConfigOption
    @Description("""
        Include the hostname of each task in the execution trace (default: `false`).
    """)
    final boolean fetchNodeName

    @ConfigOption
    @Description("""
        The FUSE device plugin to be used when enabling Fusion in unprivileged mode (default: `['nextflow.io/fuse': 1]`).
    """)
    final Map fuseDevicePlugin

    @ConfigOption
    @Description("""
        The Kubernetes HTTP client request connection timeout e.g. `'60s'`.
    """)
    final Duration httpConnectTimeout

    @ConfigOption
    @Description("""
        The Kubernetes HTTP client request connection read timeout e.g. `'60s'`.
    """)
    final Duration httpReadTimeout

    @ConfigOption
    @Description("""
        The strategy for pulling container images. Can be `IfNotPresent`, `Always`, `Never`.

        [Read more](https://kubernetes.io/docs/concepts/containers/images/#image-pull-policy)
    """)
    final String imagePullPolicy

    @ConfigOption
    @Description("""
        The path where the workflow is launched and the user data is stored (default: `<volume-claim-mount-path>/<user-name>`). Must be a path in a shared K8s persistent volume.
    """)
    final String launchDir

    @ConfigOption
    @Description("""
        The Kubernetes namespace to use (default: `default`).
    """)
    final String namespace

    @ConfigOption(types=[List, Map])
    @Description("""
        Additional pod configuration options such as environment variables, config maps, secrets, etc. Allows the same settings as the [pod](https://nextflow.io/docs/latest/process.html#pod) process directive.
    """)
    final PodOptions pod

    @ConfigOption
    @Description("""
        The path where Nextflow projects are downloaded (default: `<volume-claim-mount-path>/projects`). Must be a path in a shared K8s persistent volume.
    """)
    final String projectDir

    @Deprecated
    @ConfigOption
    @Description("""
    """)
    final String pullPolicy

    final K8sRetryConfig retryPolicy

    @ConfigOption(types=[Integer, String])
    @Description("""
        The user ID to be used to run the containers. Shortcut for the `securityContext` option.
    """)
    final Object runAsUser

    @ConfigOption
    @Description("""
        The [security context](https://kubernetes.io/docs/tasks/configure-pod-container/security-context/) to use for all pods.
    """)
    final Map securityContext

    @ConfigOption
    @Description("""
        The Kubernetes [service account name](https://kubernetes.io/docs/tasks/configure-pod-container/configure-service-account/) to use.
    """)
    final String serviceAccount

    @ConfigOption
    @Description("""
        The name of the persistent volume claim where the shared work directory is stored.
    """)
    final String storageClaimName

    @ConfigOption
    @Description("""
        The mount path for the persistent volume claim (default: `/workspace`).
    """)
    final String storageMountPath

    @ConfigOption
    @Description("""
        The path in the persistent volume to be mounted (default: `/`).
    """)
    final String storageSubPath

    @ConfigOption
    @Description("""
    """)
    final String userName

    @ConfigOption
    @Description("""
        The path of the shared work directory (default: `<user-dir>/work`). Must be a path in a shared K8s persistent volume.
    """)
    final String workDir

    /* required by extension point -- do not remove */
    K8sConfig() {
        this(Collections.emptyMap())
    }

    K8sConfig(Map opts) {
        autoMountHostPaths = opts.autoMountHostPaths as boolean
        cleanup = opts.cleanup as Boolean
        client = opts.client as Map
        clientRefreshInterval = opts.clientRefreshInterval as Duration ?: Duration.of('50m')
        clientCache = CacheBuilder.newBuilder()
            .expireAfterWrite(clientRefreshInterval.toMillis(), TimeUnit.MILLISECONDS)
            .build()
        computeResourceType = opts.computeResourceType
        context = opts.context
        cpuLimits = opts.cpuLimits as boolean
        debug = new K8sDebug(opts.debug as Map ?: Collections.emptyMap())
        fetchNodeName = opts.fetchNodeName as boolean
        fuseDevicePlugin = parseFuseDevicePlugin(opts.fuseDevicePlugin)
        httpConnectTimeout = opts.httpConnectTimeout as Duration
        httpReadTimeout = opts.httpReadTimeout as Duration
        imagePullPolicy = opts.pullPolicy ?: opts.imagePullPolicy
        namespace = opts.namespace
        pod = createPodOptions(opts.pod)
        retryPolicy = new K8sRetryConfig(opts.retryPolicy as Map ?: Collections.emptyMap())
        runAsUser = opts.runAsUser
        securityContext = opts.securityContext as Map
        serviceAccount = opts.serviceAccount
        storageClaimName = opts.storageClaimName
        storageMountPath = opts.storageMountPath ?: '/workspace'
        storageSubPath = opts.storageSubPath
        userName = opts.userName

        launchDir = opts.launchDir ?: "${storageMountPath}/${getUserName()}"
        projectDir = opts.projectDir ?: "${storageMountPath}/projects"
        workDir = opts.workDir ?: "${launchDir}/work"

        // -- shortcut to pod image pull-policy
        if( imagePullPolicy )
            pod.imagePullPolicy = imagePullPolicy

        // -- shortcut to pod volume claim
        if( storageClaimName ) {
            final volumeClaim = new PodVolumeClaim(storageClaimName, storageMountPath, storageSubPath)
            pod.volumeClaims.add(volumeClaim)
        }

        // -- shortcut to pod security context
        if( runAsUser )
            pod.securityContext = new PodSecurityContext(runAsUser)
        else if( securityContext )
            pod.securityContext = new PodSecurityContext(securityContext)
    }

    private PodOptions createPodOptions( value ) {
        if( value instanceof List )
            return new PodOptions( value as List )

        if( value instanceof Map )
            return new PodOptions( [(Map)value] )

        if( value == null )
            return new PodOptions()

        throw new IllegalArgumentException("Not a valid pod setting: $value")
    }

    Map<String,String> getLabels() {
        pod.getLabels()
    }

    Map<String,String> getAnnotations() {
        pod.getAnnotations()
    }

    boolean getCleanup(boolean defValue=true) {
        cleanup == null ? defValue : cleanup
    }

    String getUserName() {
        userName ?: System.properties.get('user.name')
    }

    Map<String,?> fuseDevicePlugin() {
        fuseDevicePlugin
    }

    Map<String,?> parseFuseDevicePlugin(Object value) {
        if( value instanceof Map && value.size()==1 )
            return value as Map<String,?>
        if( value != null )
            log.warn1 "Setting 'k8s.fuseDevicePlugin' should be a map containing exactly one entry - offending value: $value"
        return DEFAULT_FUSE_PLUGIN
    }

    /**
     * Whenever the pod should honour the entrypoint defined by the image (default: false)
     *
     *  @return  When {@code false} the launcher script is run by using pod `command` attributes which
     *      overrides the entrypoint point defined by the image.
     *
     *      When {@code true} the launcher is run via the pod `args` attribute, without altering the
     *      container entrypoint (it does however require to have a bash shell as the image entrypoint)
     *
     */
    boolean entrypointOverride() {
        return ContainerHelper.entrypointOverride()
    }

    boolean useJobResource() { ResourceType.Job.name() == computeResourceType }

    String getNextflowImageName() {
        return "nextflow/nextflow:${BuildInfo.version}"
    }

    PodOptions getPodOptions() {
        pod
    }

    boolean fetchNodeName() {
        fetchNodeName
    }

    /**
     * @return the collection of defined volume claim names
     */
    Collection<String> getClaimNames() {
        pod.volumeClaims.collect { it.claimName }
    }

    Collection<String> getClaimPaths() {
        pod.volumeClaims.collect { it.mountPath }
    }

    boolean cpuLimitsEnabled() {
        cpuLimits
    }

    /**
     * Find a volume claim name given the mount path
     *
     * @param path The volume claim mount path
     * @return The volume claim name for the given mount path
     */
    String findVolumeClaimByPath(String path) {
        final result = pod.volumeClaims.find { path.startsWith(it.mountPath) }
        return result ? result.claimName : null
    }

    ClientConfig getClient() {
        return clientCache.get('client', this::getClient0)
    }

    private ClientConfig getClient0() {
        final result = client != null
                ? clientFromNextflow(client, namespace, serviceAccount)
                : clientDiscovery(context, namespace, serviceAccount)

        if( httpConnectTimeout )
            result.httpConnectTimeout = httpConnectTimeout

        if( httpReadTimeout )
            result.httpReadTimeout = httpReadTimeout

        if( retryPolicy )
            result.retryConfig = retryPolicy

        return result
    }

    /**
     * Get the K8s client config from the declaration made in the Nextflow config file
     *
     * @param map
     *      A map representing the clint configuration options define in the nextflow
     *      config file
     * @param namespace
     *      The K8s namespace to be used. If omitted {@code default} is used.
     * @param serviceAccount
     *      The K8s service account to be used. If omitted {@code default} is used.
     * @return
     *      The Kubernetes {@link ClientConfig} object
     */
    @PackageScope ClientConfig clientFromNextflow(Map map, @Nullable String namespace, @Nullable String serviceAccount ) {
        ClientConfig.fromNextflowConfig(map,namespace,serviceAccount)
    }

    /**
     * Discover the K8s client config from the execution environment
     * that can be either a `.kube/config` file or service meta file
     * when running in a pod.
     *
     * @param contextName
     *      The name of the configuration context to be used
     * @param namespace
     *      The Kubernetes namespace to be used
     * @param serviceAccount
     *      The Kubernetes serviceAccount to be used
     * @return
     *      The discovered Kube {@link ClientConfig} object
     */
    @PackageScope ClientConfig clientDiscovery(String contextName, String namespace, String serviceAccount) {
        ClientConfig.discover(contextName, namespace, serviceAccount)
    }

    void checkStorageAndPaths(K8sClient client) {
        if( !storageClaimName )
            throw new AbortOperationException("Missing K8s storage volume claim -- The name of a persistence volume claim needs to be provided in the nextflow configuration file")

        log.debug "Kubernetes workDir=$workDir; projectDir=$projectDir; volumeClaims=${getClaimNames()}"

        for( String name : getClaimNames() ) {
            try {
                client.volumeClaimRead(name)
            }
            catch (K8sResponseException e) {
                if( e.response.code == 404 ) {
                    throw new AbortOperationException("Unknown volume claim: $name -- make sure a persistent volume claim with the specified name is defined in your K8s cluster")
                }
                else throw e
            }
        }

        if( !findVolumeClaimByPath(launchDir) )
            throw new AbortOperationException("Kubernetes `launchDir` must be a path mounted as a persistent volume -- launchDir=$launchDir; volumes=${getClaimPaths().join(', ')}")

        if( !findVolumeClaimByPath(workDir) )
            throw new AbortOperationException("Kubernetes `workDir` must be a path mounted as a persistent volume -- workDir=$workDir; volumes=${getClaimPaths().join(', ')}")

        if( !findVolumeClaimByPath(projectDir) )
            throw new AbortOperationException("Kubernetes `projectDir` must be a path mounted as a persistent volume -- projectDir=$projectDir; volumes=${getClaimPaths().join(', ')}")

    }

    static class K8sDebug implements ConfigScope {

        @ConfigOption
        @Description("""
            Save the pod spec for each task to `.command.yaml` in the task directory (default: `false`).
        """)
        final boolean yaml

        K8sDebug(Map opts) {
            yaml = opts.yaml as boolean
        }
    }
}

