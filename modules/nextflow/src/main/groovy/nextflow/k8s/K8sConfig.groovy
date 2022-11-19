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

package nextflow.k8s

import javax.annotation.Nullable

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Const
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
@Slf4j
@CompileStatic
class K8sConfig implements Map<String,Object> {

    @Delegate
    private Map<String,Object> target

    private PodOptions podOptions

    K8sConfig(Map<String,Object> config) {
        target = config ?: Collections.<String,Object>emptyMap()

        this.podOptions = createPodOptions(target.pod)
        if( getStorageClaimName() ) {
            final name = getStorageClaimName()
            final mount = getStorageMountPath()
            final subPath = getStorageSubPath()
            this.podOptions.volumeClaims.add(new PodVolumeClaim(name, mount, subPath))
        }

        // -- shortcut to pod image pull-policy
        if( target.pullPolicy )
            podOptions.imagePullPolicy = target.pullPolicy.toString()
        else if( target.imagePullPolicy )
            podOptions.imagePullPolicy = target.imagePullPolicy.toString()

        // -- shortcut to pod security context
        if( target.runAsUser != null )
            podOptions.securityContext = new PodSecurityContext(target.runAsUser)
        else if( target.securityContext instanceof Map )
            podOptions.securityContext = new PodSecurityContext(target.securityContext as Map)
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
        podOptions.getLabels()
    }

    Map<String,String> getAnnotations() {
        podOptions.getAnnotations()
    }

    K8sDebug getDebug() {
        new K8sDebug( (Map<String,Object>)get('debug') )
    }

    boolean getCleanup(boolean defValue=true) {
        target.cleanup == null ? defValue : Boolean.valueOf( target.cleanup as String )
    }

    String getUserName() {
        target.userName ?: System.properties.get('user.name')
    }

    String getStorageClaimName() {
        target.storageClaimName as String
    }

    String getStorageMountPath() {
        target.storageMountPath ?: '/workspace' as String
    }

    String getStorageSubPath() {
        target.storageSubPath
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
        def result = target.entrypointOverride
        if( result == null )
            result = System.getenv('NXF_CONTAINER_ENTRYPOINT_OVERRIDE')
        return result
    }

    /**
     * @return the path where the workflow is launched and the user data is stored
     */
    String getLaunchDir() {
        if( target.userDir ) {
            log.warn "K8s `userDir` has been deprecated -- Use `launchDir` instead"
            return target.userDir
        }
        target.launchDir ?: "${getStorageMountPath()}/${getUserName()}" as String
    }

    /**
     * @return Defines the path where the workflow temporary data is stored. This must be a path in a shared K8s persistent volume (default:<user-dir>/work).
     */
    String getWorkDir() {
        target.workDir ?: "${getLaunchDir()}/work" as String
    }

    /**
     * @return Defines the path where Nextflow projects are downloaded. This must be a path in a shared K8s persistent volume (default: <volume-claim-mount-path>/projects).
     */
    String getProjectDir() {
        target.projectDir ?: "${getStorageMountPath()}/projects" as String
    }

    String getNamespace() { target.namespace }

    boolean useJobResource() { ResourceType.Job.name() == target.computeResourceType?.toString() }

    String getServiceAccount() { target.serviceAccount }

    String getNextflowImageName() {
        final defImage = "nextflow/nextflow:${Const.APP_VER}"
        return target.navigate('nextflow.image', defImage)
    }

    boolean getAutoMountHostPaths() {
        Boolean.valueOf( target.autoMountHostPaths as String )
    }

    PodOptions getPodOptions() {
        podOptions
    }

    @Memoized
    boolean fetchNodeName() {
        Boolean.valueOf( target.fetchNodeName as String )
    }

    /**
     * @return the collection of defined volume claim names
     */
    Collection<String> getClaimNames() {
        podOptions.volumeClaims.collect { it.claimName }
    }

    Collection<String> getClaimPaths() {
        podOptions.volumeClaims.collect { it.mountPath }
    }

    /**
     * Find a volume claim name given the mount path
     *
     * @param path The volume claim mount path
     * @return The volume claim name for the given mount path
     */
    String findVolumeClaimByPath(String path) {
        def result = podOptions.volumeClaims.find { path.startsWith(it.mountPath) }
        return result ? result.claimName : null
    }

    @Memoized
    ClientConfig getClient() {

        final result = ( target.client instanceof Map
                ? clientFromNextflow(target.client as Map, target.namespace as String, target.serviceAccount as String)
                : clientDiscovery(target.context as String, target.namespace as String, target.serviceAccount as String)
        )

        if( target.httpConnectTimeout )
            result.httpConnectTimeout = target.httpConnectTimeout as Duration

        if( target.httpReadTimeout )
            result.httpReadTimeout = target.httpReadTimeout as Duration

        if( target.maxErrorRetry )
            result.maxErrorRetry = target.maxErrorRetry as Integer

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
        if( !getStorageClaimName() )
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

        if( !findVolumeClaimByPath(getLaunchDir()) )
            throw new AbortOperationException("Kubernetes `launchDir` must be a path mounted as a persistent volume -- launchDir=$launchDir; volumes=${getClaimPaths().join(', ')}")

        if( !findVolumeClaimByPath(getWorkDir()) )
            throw new AbortOperationException("Kubernetes `workDir` must be a path mounted as a persistent volume -- workDir=$workDir; volumes=${getClaimPaths().join(', ')}")

        if( !findVolumeClaimByPath(getProjectDir()) )
            throw new AbortOperationException("Kubernetes `projectDir` must be a path mounted as a persistent volume -- projectDir=$projectDir; volumes=${getClaimPaths().join(', ')}")


    }

    @CompileStatic
    static class K8sDebug {

        @Delegate
        Map<String,Object> target

        K8sDebug(Map<String,Object> debug) {
            this.target = debug ?: Collections.<String,Object>emptyMap()
        }

        boolean getYaml() { Boolean.valueOf( target.yaml as String ) }
    }
}

