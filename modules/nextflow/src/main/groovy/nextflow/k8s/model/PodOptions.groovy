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

package nextflow.k8s.model

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope
import groovy.transform.ToString

/**
 * Model K8s pod options such as environment variables,
 * secret and config-maps
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true)
@EqualsAndHashCode(includeFields = true)
class PodOptions {

    private String imagePullPolicy

    private String imagePullSecret

    private Collection<PodEnv> envVars

    private Collection<PodMountConfig> mountConfigMaps

    private Collection<PodMountCsiEphemeral> mountCsiEphemerals

    private Collection<PodMountEmptyDir> mountEmptyDirs

    private Collection<PodMountSecret> mountSecrets

    private Collection<PodVolumeClaim> mountClaims

    private Collection<PodHostMount> mountHostPaths

    private Map<String,String> labels = [:]

    private Map<String,String> annotations = [:]

    private PodNodeSelector nodeSelector

    private Map affinity

    private PodSecurityContext securityContext

    private boolean automountServiceAccountToken

    private String priorityClassName

    private List<Map> tolerations

    private Boolean privileged

    private String schedulerName

    private Integer ttlSecondsAfterFinished
    
    PodOptions( List<Map> options=null ) {
        int size = options ? options.size() : 0
        envVars = new HashSet<>(size)
        mountConfigMaps = new HashSet<>(size)
        mountCsiEphemerals = new HashSet<>(size)
        mountEmptyDirs = new HashSet<>(size)
        mountSecrets = new HashSet<>(size)
        mountClaims = new HashSet<>(size)
        mountHostPaths = new HashSet<>(10)
        automountServiceAccountToken = true
        tolerations = new ArrayList<Map>(size)
        init(options)
    }

    @PackageScope void init(List<Map> options) {
        if( !options ) return 
        for( Map entry : options ) {
            create(entry)
        }
    }

    @PackageScope void create(Map<String,String> entry) {
        if( entry.env && entry.value ) {
            envVars << PodEnv.value(entry.env, entry.value)
        }
        else if( entry.env && entry.fieldPath ) {
            envVars << PodEnv.fieldPath(entry.env, entry.fieldPath)
        }
        else if( entry.env && entry.secret ) {
            envVars << PodEnv.secret(entry.env, entry.secret)
        }
        else if( entry.env && entry.config ) {
            envVars << PodEnv.config(entry.env, entry.config)
        }
        else if( entry.mountPath && entry.secret ) {
            mountSecrets << new PodMountSecret(entry)
        }
        else if( entry.mountPath && entry.config ) {
            mountConfigMaps << new PodMountConfig(entry)
        }
        else if( entry.mountPath && entry.csi ) {
            mountCsiEphemerals << new PodMountCsiEphemeral(entry)
        }
        else if( entry.mountPath && entry.emptyDir != null ) {
            mountEmptyDirs << new PodMountEmptyDir(entry)
        }
        else if( entry.mountPath && entry.volumeClaim ) {
            mountClaims << new PodVolumeClaim(entry)
        }
        else if( entry.mountPath && entry.hostPath instanceof CharSequence ) {
            mountHostPaths << new PodHostMount(entry.hostPath, entry.mountPath)
        }
        else if( entry.pullPolicy || entry.imagePullPolicy ) {
            this.imagePullPolicy = entry.pullPolicy ?: entry.imagePullPolicy as String
        }
        else if( entry.imagePullSecret || entry.imagePullSecrets ) {
            this.imagePullSecret = entry.imagePullSecret ?: entry.imagePullSecrets
        }
        else if( entry.label && entry.value ) {
            this.labels.put(entry.label as String, entry.value as String)
        }
        else if( entry.runAsUser != null ) {
            this.securityContext = new PodSecurityContext(entry.runAsUser)
        }
        else if( entry.securityContext instanceof Map ) {
            this.securityContext = new PodSecurityContext(entry.securityContext as Map)
        }
        else if( entry.nodeSelector ) {
            this.nodeSelector = new PodNodeSelector(entry.nodeSelector)
        }
        else if( entry.affinity instanceof Map ) {
            this.affinity = entry.affinity as Map
        }
        else if( entry.annotation && entry.value ) {
            this.annotations.put(entry.annotation as String, entry.value as String)
        }
        else if( entry.automountServiceAccountToken instanceof Boolean ) {
            this.automountServiceAccountToken = entry.automountServiceAccountToken as Boolean
        }
        else if( entry.priorityClassName ) {
            this.priorityClassName = entry.priorityClassName
        }
        else if( entry.toleration instanceof Map ) {
            tolerations << (entry.toleration as Map)
        }
        else if( entry.privileged instanceof Boolean ) {
            this.privileged = entry.privileged as Boolean
        }
        else if( entry.schedulerName ) {
            this.schedulerName = entry.schedulerName
        }
        else if( entry.ttlSecondsAfterFinished instanceof Integer ) {
            this.ttlSecondsAfterFinished = entry.ttlSecondsAfterFinished as Integer
        }
        else
            throw new IllegalArgumentException("Unknown pod options: $entry")
    }


    Collection<PodEnv> getEnvVars() { envVars }

    Collection<PodMountConfig> getMountConfigMaps() { mountConfigMaps }

    Collection<PodMountCsiEphemeral> getMountCsiEphemerals() { mountCsiEphemerals }

    Collection<PodMountEmptyDir> getMountEmptyDirs() { mountEmptyDirs }

    Collection<PodMountSecret> getMountSecrets() { mountSecrets }

    Collection<PodHostMount> getMountHostPaths() { mountHostPaths }

    Collection<PodVolumeClaim> getVolumeClaims() { mountClaims }

    Map<String,String> getLabels() { labels }

    Map<String,String> getAnnotations() { annotations }

    PodNodeSelector getNodeSelector() { nodeSelector }

    PodOptions setNodeSelector( PodNodeSelector sel ) {
        nodeSelector = sel
        return this
    }

    Map getAffinity() { affinity }

    PodSecurityContext getSecurityContext() { securityContext }

    PodOptions setSecurityContext( PodSecurityContext ctx ) {
        this.securityContext = ctx
        return this
    }

    String getImagePullSecret() { imagePullSecret }

    PodOptions setImagePullSecret( String secret ) {
        this.imagePullSecret = secret
        return this
    }

    String getImagePullPolicy() { imagePullPolicy }

    PodOptions setImagePullPolicy( String policy ) {
        this.imagePullPolicy = policy
        return this
    }

    boolean getAutomountServiceAccountToken() { automountServiceAccountToken }

    PodOptions setAutomountServiceAccountToken( boolean mount ) {
        this.automountServiceAccountToken = mount
        return this
    }

    String getPriorityClassName() { priorityClassName }

    String getSchedulerName() { schedulerName }

    List<Map> getTolerations() { tolerations }

    Boolean getPrivileged() { privileged }

    Integer getTtlSecondsAfterFinished() { ttlSecondsAfterFinished }

    PodOptions plus( PodOptions other ) {
        def result = new PodOptions()

        // env vars
        result.envVars.addAll(envVars)
        result.envVars.addAll( other.envVars )

        // config maps
        result.mountConfigMaps.addAll( mountConfigMaps )
        result.mountConfigMaps.addAll( other.mountConfigMaps )

        // csi ephemeral volumes
        result.mountCsiEphemerals.addAll( mountCsiEphemerals )
        result.mountCsiEphemerals.addAll( other.mountCsiEphemerals )

        // empty dirs
        result.mountEmptyDirs.addAll( mountEmptyDirs )
        result.mountEmptyDirs.addAll( other.mountEmptyDirs )

        // host paths
        result.mountHostPaths.addAll( mountHostPaths )
        result.mountHostPaths.addAll( other.mountHostPaths )

        // secrets
        result.mountSecrets.addAll( mountSecrets )
        result.mountSecrets.addAll( other.mountSecrets )

        // volume claims
        result.volumeClaims.addAll( volumeClaims )
        result.volumeClaims.addAll( other.volumeClaims )

        // sec context
        result.securityContext = other.securityContext ?: this.securityContext

        // node selector
        result.nodeSelector = other.nodeSelector ?: this.nodeSelector

        // affinity
        result.affinity = other.affinity ?: this.affinity

        // pull policy
        result.imagePullPolicy = other.imagePullPolicy ?: this.imagePullPolicy

        // image secret
        result.imagePullSecret = other.imagePullSecret ?: this.imagePullSecret

        // labels
        result.labels.putAll(labels)
        result.labels.putAll(other.labels)

        // annotations
        result.annotations.putAll(annotations)
        result.annotations.putAll(other.annotations)

        // automount service account token
        result.automountServiceAccountToken = other.automountServiceAccountToken & this.automountServiceAccountToken

        // priority class name
        result.priorityClassName = other.priorityClassName ?: this.priorityClassName

        // tolerations
        result.tolerations = other.tolerations ?: this.tolerations

        //  privileged execution
        result.privileged = other.privileged!=null ? other.privileged : this.privileged

        // scheduler name
        result.schedulerName = other.schedulerName ?: this.schedulerName

        // ttl seconds after finished (job)
        result.ttlSecondsAfterFinished = other.ttlSecondsAfterFinished!=null ? other.ttlSecondsAfterFinished : this.ttlSecondsAfterFinished

        return result
    }
}
