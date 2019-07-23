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

    private Collection<PodMountSecret> mountSecrets

    private Collection<PodVolumeClaim> mountClaims

    private Map<String,String> labels = [:]

    private Map<String,String> annotations = [:]

    private PodNodeSelector nodeSelector

    private PodSecurityContext securityContext

    PodOptions( List<Map> options=null ) {
        int size = options ? options.size() : 0
        envVars = new HashSet<>(size)
        mountSecrets = new HashSet<>(size)
        mountConfigMaps = new HashSet<>(size)
        mountClaims = new HashSet<>(size)
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
        else if( entry.env && entry.secret ) {
            envVars << PodEnv.secret(entry.env, entry.secret)
        }
        else if( entry.env && entry.config ) {
            envVars << PodEnv.config(entry.env, entry.config)
        }
        else if( entry.mountPath && entry.secret ) {
            mountSecrets <<  new PodMountSecret(entry)
        }
        else if( entry.mountPath && entry.config ) {
            mountConfigMaps << new PodMountConfig(entry)
        }
        else if( entry.mountPath && entry.volumeClaim ) {
            mountClaims << new PodVolumeClaim(entry)
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
        else if( entry.annotation && entry.value ) {
            this.annotations.put(entry.annotation as String, entry.value as String)
        }
        else 
            throw new IllegalArgumentException("Unknown pod options: $entry")
    }


    Collection<PodEnv> getEnvVars() { envVars }

    Collection<PodMountConfig> getMountConfigMaps() { mountConfigMaps }

    Collection<PodMountSecret> getMountSecrets() { mountSecrets }

    Collection<PodVolumeClaim> getVolumeClaims() { mountClaims }

    Map<String,String> getLabels() { labels }

    Map<String,String> getAnnotations() { annotations }

    PodSecurityContext getSecurityContext() { securityContext }

    PodNodeSelector getNodeSelector() { nodeSelector }

    PodOptions setNodeSelector( PodNodeSelector sel ) {
        nodeSelector = sel
        return this
    }

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

    PodOptions plus( PodOptions other ) {
        def result = new PodOptions()
        // env vars
        result.envVars.addAll(envVars)
        result.envVars.addAll( other.envVars )

        // config maps
        result.mountConfigMaps.addAll( mountConfigMaps )
        result.mountConfigMaps.addAll( other.mountConfigMaps )

        // secrets
        result.mountSecrets.addAll( mountSecrets )
        result.mountSecrets.addAll( other.mountSecrets )

        // volume claims
        result.volumeClaims.addAll( volumeClaims )
        result.volumeClaims.addAll( other.volumeClaims )

        // sec context
        if( other.securityContext )
            result.securityContext = other.securityContext
        else
            result.securityContext = securityContext

        // node select
        result.nodeSelector = other.nodeSelector ?: this.nodeSelector

        // pull policy
        if (other.imagePullPolicy)
            result.imagePullPolicy = other.imagePullPolicy
        else
            result.imagePullPolicy = imagePullPolicy

        // image secret
        if (other.imagePullSecret)
            result.imagePullSecret = other.imagePullSecret
        else
            result.imagePullSecret = imagePullSecret

        // labels
        result.labels.putAll(labels)
        result.labels.putAll(other.labels)

        // annotations
        result.annotations.putAll(annotations)
        result.annotations.putAll(other.annotations)

        return result
    }
}
