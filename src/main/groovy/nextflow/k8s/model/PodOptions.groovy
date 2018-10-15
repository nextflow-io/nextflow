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
        else 
            throw new IllegalArgumentException("Unknown pod options: $entry")
    }


    Collection<PodEnv> getEnvVars() { envVars }

    Collection<PodMountConfig> getMountConfigMaps() { mountConfigMaps }

    Collection<PodMountSecret> getMountSecrets() { mountSecrets }

    Collection<PodVolumeClaim> getVolumeClaims() { mountClaims }

    Map<String,String> getLabels() { labels }

    PodSecurityContext getSecurityContext() { securityContext }

    PodOptions setSecurityContext( PodSecurityContext ctx) {
        this.securityContext = ctx
        return this
    }

    String getImagePullSecret() { imagePullSecret }

    String getImagePullPolicy() { imagePullPolicy }

    PodOptions setImagePullPolicy(String p ) {
        this.imagePullPolicy = p
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

        if( other.securityContext )
            result.securityContext = other.securityContext
        else
            result.securityContext = securityContext

        return result
    }
}
