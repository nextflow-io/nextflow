/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

    K8sDebug getDebug() {
        new K8sDebug( (Map<String,Object>)get('debug') )
    }

    boolean getCleanup() {
        target.cleanup == null || target.cleanup as boolean
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

    String getServiceAccount() { target.serviceAccount }

    String getNextflowImageName() {
        final defImage = "nextflow/nextflow:${Const.APP_VER}"
        return target.navigate('nextflow.image', defImage)
    }

    boolean getAutoMountHostPaths() {
        target.autoMountHostPaths as boolean
    }

    PodOptions getPodOptions() {
        podOptions
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
                ? clientFromMap(target.client as Map)
                : clientDiscovery(target.context as String)
        )

        if( target.namespace ) {
            result.namespace = target.namespace as String
        }

        if( target.serviceAccount ) {
            result.serviceAccount = target.serviceAccount as String
        }

        return result
    }

    @PackageScope ClientConfig clientFromMap( Map map ) {
        ClientConfig.fromMap(map)
    }

    @PackageScope ClientConfig clientDiscovery( String ctx ) {
        ClientConfig.discover(ctx)
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

        boolean getYaml() { target.yaml as boolean }
    }
}

