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

import java.nio.file.Path
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.util.MemoryUnit
/**
 * Object build for a K8s pod specification
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class PodSpecBuilder {

    static @PackageScope AtomicInteger VOLUMES = new AtomicInteger()

    String podName

    String imageName

    String imagePullPolicy

    String imagePullSecret

    List<String> command = []

    Map<String,String> labels = [:]

    String namespace

    String restart

    List<PodEnv> envVars = []

    String workDir

    Integer cpus

    String memory

    String serviceAccount

    Collection<PodMountSecret> secrets = []

    Collection<PodMountConfig> configMaps = []

    Collection<PodHostMount> hostMounts = []

    Collection<PodVolumeClaim> volumeClaims = []

    PodSecurityContext securityContext

    PodNodeSelector nodeSelector

    /**
     * @return A sequential volume unique identifier
     */
    static protected String nextVolName() {
        "vol-${VOLUMES.incrementAndGet()}".toString()
    }


    PodSpecBuilder withPodName(String name) {
        this.podName = name
        return this
    }

    PodSpecBuilder withImageName(String name) {
        this.imageName = name
        return this
    }

    PodSpecBuilder withImagePullPolicy(String policy) {
        this.imagePullPolicy = policy
        return this
    }

    PodSpecBuilder withWorkDir( String path ) {
        this.workDir = path
        return this
    }

    PodSpecBuilder withWorkDir(Path path ) {
        this.workDir = path.toString()
        return this
    }

    PodSpecBuilder withNamespace(String name) {
        this.namespace = name
        return this
    }

    PodSpecBuilder withServiceAccount(String name) {
        this.serviceAccount = name
        return this
    }

    PodSpecBuilder withCommand( cmd ) {
        assert cmd instanceof List || cmd instanceof CharSequence, "Missing or invalid K8s command parameter: $cmd"
        this.command = cmd instanceof List ? cmd : ['/bin/bash','-c', cmd.toString()]
        return this
    }

    PodSpecBuilder withCpus( Integer cpus ) {
        this.cpus = cpus
        return this
    }

    PodSpecBuilder withMemory(String mem) {
        this.memory = mem
        return this
    }

    PodSpecBuilder withMemory(MemoryUnit mem)  {
        this.memory = "${mem.mega}Mi".toString()
        return this
    }

    PodSpecBuilder withLabel( String name, String value ) {
        this.labels.put(name, value)
        return this
    }

    PodSpecBuilder withLabels(Map labels) {
        this.labels.putAll(labels)
        return this
    }


    PodSpecBuilder withEnv( PodEnv var ) {
        envVars.add(var)
        return this
    }

    PodSpecBuilder withEnv( Collection vars ) {
        envVars.addAll(vars)
        return this
    }

    PodSpecBuilder withVolumeClaim( PodVolumeClaim claim ) {
        volumeClaims.add(claim)
        return this
    }

    PodSpecBuilder withVolumeClaims( Collection<PodVolumeClaim> claims ) {
        volumeClaims.addAll(claims)
        return this
    }

    PodSpecBuilder withConfigMaps( Collection<PodMountConfig> configMaps ) {
        this.configMaps.addAll(configMaps)
        return this
    }

    PodSpecBuilder withConfigMap( PodMountConfig configMap ) {
        this.configMaps.add(configMap)
        return this
    }

    PodSpecBuilder withSecrets( Collection<PodMountSecret> secrets ) {
        this.secrets.addAll(secrets)
        return this
    }

    PodSpecBuilder withSecret( PodMountSecret secret ) {
        this.secrets.add(secret)
        return this
    }

    PodSpecBuilder withHostMounts( Collection<PodHostMount> mounts ) {
        this.hostMounts.addAll(mounts)
        return this
    }

    PodSpecBuilder withHostMount( String host, String mount ) {
        this.hostMounts.add( new PodHostMount(host, mount))
        return this
    }

    PodSpecBuilder withPodOptions(PodOptions opts) {
        // -- pull policy
        if( opts.imagePullPolicy )
            imagePullPolicy = opts.imagePullPolicy
        if( opts.imagePullSecret )
            imagePullSecret = opts.imagePullSecret
        // -- env vars
        if( opts.getEnvVars() )
            envVars.addAll( opts.getEnvVars() )
        // -- secrets
        if( opts.getMountSecrets() )
            secrets.addAll( opts.getMountSecrets() )
        // -- configMaps
        if( opts.getMountConfigMaps() )
            configMaps.addAll( opts.getMountConfigMaps() )
        // -- volume claims 
        if( opts.getVolumeClaims() )
            volumeClaims.addAll( opts.getVolumeClaims() )
        // -- labels
        if( opts.labels ) {
            def keys = opts.labels.keySet()
            if( 'app' in keys ) throw new IllegalArgumentException("Invalid pod label -- `app` is a reserved label")
            if( 'runName' in keys ) throw new IllegalArgumentException("Invalid pod label -- `runName` is a reserved label")
            labels.putAll( opts.labels )
        }
        // -- security context
        if( opts.securityContext )
            securityContext = opts.securityContext
        if( opts.nodeSelector )
            nodeSelector = opts.nodeSelector

        return this
    }

    @PackageScope List<Map> createPullSecret() {
        def result = new ArrayList(1)
        def entry = new LinkedHashMap(1)
        entry.name = imagePullSecret
        result.add(entry)
        return result
    }

    Map build() {
        assert this.podName, 'Missing K8s podName parameter'
        assert this.imageName, 'Missing K8s imageName parameter'
        assert this.command, 'Missing K8s command parameter'

        final restart = this.restart ?: 'Never'

        final metadata = new LinkedHashMap<String,Object>()
        metadata.name = podName
        metadata.namespace = namespace ?: 'default'

        final labels = this.labels ?: [:]
        final env = []
        for( PodEnv entry : this.envVars ) {
            env.add(entry.toSpec())
        }

        final res = [:]
        if( this.cpus )
            res.cpu = this.cpus
        if( this.memory )
            res.memory = this.memory

        final container = [
                name: this.podName,
                image: this.imageName,
                command: this.command
        ]
        
        if( this.workDir )
            container.put('workingDir', workDir)

        if( imagePullPolicy )
            container.imagePullPolicy = imagePullPolicy

        final spec = [
                restartPolicy: restart,
                containers: [ container ],
        ]

        if( nodeSelector )
            spec.nodeSelector = nodeSelector.toSpec()

        if( this.serviceAccount )
            spec.serviceAccountName = this.serviceAccount

        if( securityContext )
            spec.securityContext = securityContext.toSpec()

        if( imagePullSecret )
            spec.imagePullSecrets = createPullSecret()

        // add labels
        if( labels )
            metadata.labels = labels

        final pod = [
                apiVersion: 'v1',
                kind: 'Pod',
                metadata: metadata,
                spec: spec
        ]

        // add environment
        if( env )
            container.env = env

        // add resources
        if( res ) {
            Map limits = [limits: res]
            container.resources = limits
        }

        // add storage definitions ie. volumes and mounts
        final mounts = []
        final volumes = []

        // -- volume claims
        for( PodVolumeClaim entry : volumeClaims ) {
            final name = nextVolName()
            final claim = [name: name, mountPath: entry.mountPath ]
            if( entry.subPath ) claim.subPath = entry.subPath
            mounts << claim
            volumes << [name: name, persistentVolumeClaim: [claimName: entry.claimName]]
        }

        // -- configMap volumes
        for( PodMountConfig entry : configMaps ) {
            final name = nextVolName()
            configMapToSpec(name, entry, mounts, volumes)
        }

        // host mounts
        for( PodHostMount entry : hostMounts ) {
            final name = nextVolName()
            mounts << [name: name, mountPath: entry.mountPath]
            volumes << [name: name, hostPath: [path: entry.hostPath]]
        }

        // secret volumes
        for( PodMountSecret entry : secrets ) {
            final name = nextVolName()
            secretToSpec(name, entry, mounts, volumes)
        }


        if( volumes )
            spec.volumes = volumes
        if( mounts )
            container.volumeMounts = mounts

        return pod
    }

    @PackageScope
    @CompileDynamic
    static void secretToSpec(String volName, PodMountSecret entry, List mounts, List volumes ) {
        assert entry

        final secret = [secretName: entry.secretName]
        if( entry.secretKey ) {
            secret.items = [ [key: entry.secretKey, path: entry.fileName ] ]
        }

        mounts << [name: volName, mountPath: entry.mountPath]
        volumes << [name: volName, secret: secret ]
    }

    @PackageScope
    @CompileDynamic
    static void configMapToSpec(String volName, PodMountConfig entry, List<Map> mounts, List<Map> volumes ) {
        assert entry

        final config = [name: entry.configName]
        if( entry.configKey ) {
            config.items = [ [key: entry.configKey, path: entry.fileName ] ]
        }

        mounts << [name: volName, mountPath: entry.mountPath]
        volumes << [name: volName, configMap: config ]
    }


}
