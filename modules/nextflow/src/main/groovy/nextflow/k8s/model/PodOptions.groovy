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

    private Map affinity

    private Map<String,String> annotations = [:]

    private boolean automountServiceAccountToken

    private Collection<PodEnv> envVars

    private String imagePullPolicy

    private String imagePullSecret

    private Map<String,String> labels = [:]

    private Collection<PodMountConfig> mountConfigMaps

    private Collection<PodMountSecret> mountSecrets

    private Collection<PodVolumeClaim> mountVolumeClaims

    private PodNodeSelector nodeSelector

    private String priorityClassName

    private PodSecurityContext securityContext

    PodOptions( List<Map> options=null ) {
        int size = options ? options.size() : 0
        automountServiceAccountToken = true
        envVars = new HashSet<>(size)
        mountConfigMaps = new HashSet<>(size)
        mountSecrets = new HashSet<>(size)
        mountVolumeClaims = new HashSet<>(size)
        init(options)
    }

    @PackageScope void init(List<Map> options) {
        if( !options ) return 
        for( Map entry : options ) {
            create(entry)
        }
    }

    @PackageScope void create(Map<String,String> entry) {
        if( entry.affinity instanceof Map ) {
            this.affinity = entry.affinity as Map
        }
        else if( entry.annotation && entry.value ) {
            this.annotations.put(entry.annotation as String, entry.value as String)
        }
        else if( entry.automountServiceAccountToken instanceof Boolean ) {
            this.automountServiceAccountToken = entry.automountServiceAccountToken as Boolean
        }
        else if( entry.env && entry.config ) {
            envVars << PodEnv.config(entry.env, entry.config)
        }
        else if( entry.env && entry.fieldPath ) {
            envVars << PodEnv.fieldPath(entry.env, entry.fieldPath)
        }
        else if( entry.env && entry.secret ) {
            envVars << PodEnv.secret(entry.env, entry.secret)
        }
        else if( entry.env && entry.value ) {
            envVars << PodEnv.value(entry.env, entry.value)
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
        else if( entry.mountPath && entry.config ) {
            mountConfigMaps << new PodMountConfig(entry)
        }
        else if( entry.mountPath && entry.secret ) {
            mountSecrets << new PodMountSecret(entry)
        }
        else if( entry.mountPath && entry.volumeClaim ) {
            mountVolumeClaims << new PodVolumeClaim(entry)
        }
        else if( entry.nodeSelector ) {
            this.nodeSelector = new PodNodeSelector(entry.nodeSelector)
        }
        else if( entry.priorityClassName ) {
            this.priorityClassName = entry.priorityClassName
        }
        else if( entry.runAsUser != null ) {
            this.securityContext = new PodSecurityContext(entry.runAsUser)
        }
        else if( entry.securityContext instanceof Map ) {
            this.securityContext = new PodSecurityContext(entry.securityContext as Map)
        }
        else 
            throw new IllegalArgumentException("Unknown pod options: $entry")
    }


    Map getAffinity() { affinity }

    Map<String,String> getAnnotations() { annotations }

    boolean getAutomountServiceAccountToken() { automountServiceAccountToken }

    PodOptions setAutomountServiceAccountToken( boolean mount ) {
        this.automountServiceAccountToken = mount
        return this
    }

    Collection<PodEnv> getEnvVars() { envVars }

    String getImagePullPolicy() { imagePullPolicy }

    PodOptions setImagePullPolicy( String policy ) {
        this.imagePullPolicy = policy
        return this
    }

    String getImagePullSecret() { imagePullSecret }

    PodOptions setImagePullSecret( String secret ) {
        this.imagePullSecret = secret
        return this
    }

    Map<String,String> getLabels() { labels }

    Collection<PodMountConfig> getMountConfigMaps() { mountConfigMaps }

    Collection<PodMountSecret> getMountSecrets() { mountSecrets }

    Collection<PodVolumeClaim> getVolumeClaims() { mountVolumeClaims }

    PodNodeSelector getNodeSelector() { nodeSelector }

    PodOptions setNodeSelector( PodNodeSelector sel ) {
        nodeSelector = sel
        return this
    }

    String getPriorityClassName() { priorityClassName }

    PodSecurityContext getSecurityContext() { securityContext }

    PodOptions setSecurityContext( PodSecurityContext ctx ) {
        this.securityContext = ctx
        return this
    }

    PodOptions plus( PodOptions other ) {
        def result = new PodOptions()

        // affinity
        result.affinity = other.affinity ?: this.affinity

        // annotations
        result.annotations.putAll(annotations)
        result.annotations.putAll(other.annotations)

        // automount service account token
        result.automountServiceAccountToken = other.automountServiceAccountToken & this.automountServiceAccountToken

        // environment variables
        result.envVars.addAll(envVars)
        result.envVars.addAll( other.envVars )

        // image pull policy
        if (other.imagePullPolicy)
            result.imagePullPolicy = other.imagePullPolicy
        else
            result.imagePullPolicy = imagePullPolicy

        // image pull secret
        if (other.imagePullSecret)
            result.imagePullSecret = other.imagePullSecret
        else
            result.imagePullSecret = imagePullSecret

        // labels
        result.labels.putAll(labels)
        result.labels.putAll(other.labels)

        // config maps
        result.mountConfigMaps.addAll( mountConfigMaps )
        result.mountConfigMaps.addAll( other.mountConfigMaps )

        // secrets
        result.mountSecrets.addAll( mountSecrets )
        result.mountSecrets.addAll( other.mountSecrets )

        // volume claims
        result.volumeClaims.addAll( volumeClaims )
        result.volumeClaims.addAll( other.volumeClaims )

        // node selector
        result.nodeSelector = other.nodeSelector ?: this.nodeSelector

        // priority class name
        result.priorityClassName = other.priorityClassName ?: this.priorityClassName

        // security context
        if( other.securityContext )
            result.securityContext = other.securityContext
        else
            result.securityContext = securityContext

        return result
    }
}
