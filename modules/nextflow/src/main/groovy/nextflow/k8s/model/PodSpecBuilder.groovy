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

import java.nio.file.Path
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.executor.res.AcceleratorResource
import nextflow.util.MemoryUnit
import groovy.util.logging.Slf4j

/**
 * Object build for a K8s pod specification
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Slf4j
class PodSpecBuilder {

    static @PackageScope AtomicInteger VOLUMES = new AtomicInteger()

    String podName

    AcceleratorResource accelerator

    Map affinity

    Map<String,String> annotations = [:]

    boolean automountServiceAccountToken = true

    List<String> command = []

    Integer cpus

    List<PodEnv> envVars = []

    String imageName

    String imagePullPolicy

    String imagePullSecret

    Map<String,String> labels = [:]

    String memory

    Collection<PodMountSecret> secrets = []

    Collection<PodMountConfig> configMaps = []

    Collection<PodMountHostPath> hostPaths = []

    Collection<PodVolumeClaim> volumeClaims = []

    String namespace

    PodNodeSelector nodeSelector

    String priorityClassName

    String restartPolicy

    PodSecurityContext securityContext

    String serviceAccount

    String workDir

    /**
     * @return A sequential volume unique identifier
     */
    static protected String nextVolumeName() {
        "vol-${VOLUMES.incrementAndGet()}".toString()
    }


    PodSpecBuilder withPodName(String name) {
        this.podName = name
        return this
    }

    PodSpecBuilder withAccelerator(AcceleratorResource acc) {
        this.accelerator = acc
        return this
    }

    PodSpecBuilder withAnnotation( String name, String value ) {
        this.annotations.put(name, value)
        return this
    }

    PodSpecBuilder withAnnotations(Map annotations) {
        this.annotations.putAll(annotations)
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

    PodSpecBuilder withEnv( PodEnv env ) {
        envVars.add(env)
        return this
    }

    PodSpecBuilder withEnv( Collection envs ) {
        envVars.addAll(envs)
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

    PodSpecBuilder withLabel( String name, String value ) {
        this.labels.put(name, value)
        return this
    }

    PodSpecBuilder withLabels(Map labels) {
        this.labels.putAll(labels)
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

    PodSpecBuilder withConfigMaps( Collection<PodMountConfig> configMaps ) {
        this.configMaps.addAll(configMaps)
        return this
    }

    PodSpecBuilder withConfigMap( PodMountConfig configMap ) {
        this.configMaps.add(configMap)
        return this
    }

    PodSpecBuilder withHostPaths( Collection<PodMountHostPath> hostPaths ) {
        this.hostPaths.addAll(hostPaths)
        return this
    }

    PodSpecBuilder withHostPath( String host, String mount ) {
        this.hostPaths.add( new PodMountHostPath(host, mount))
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

    PodSpecBuilder withVolumeClaim( PodVolumeClaim claim ) {
        volumeClaims.add(claim)
        return this
    }

    PodSpecBuilder withVolumeClaims( Collection<PodVolumeClaim> claims ) {
        volumeClaims.addAll(claims)
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

    PodSpecBuilder withWorkDir( String path ) {
        this.workDir = path
        return this
    }

    PodSpecBuilder withWorkDir(Path path ) {
        this.workDir = path.toString()
        return this
    }

    PodSpecBuilder withPodOptions(PodOptions opts) {
        // -- affinity
        if( opts.affinity )
            affinity = opts.affinity

        // - annotations
        if( opts.annotations ) {
            annotations.putAll( opts.annotations )
        }

        // -- automount service account token
        automountServiceAccountToken = opts.automountServiceAccountToken

        // -- environment variables
        if( opts.getEnvVars() )
            envVars.addAll( opts.getEnvVars() )

        // -- image pull policy
        if( opts.imagePullPolicy )
            imagePullPolicy = opts.imagePullPolicy
        if( opts.imagePullSecret )
            imagePullSecret = opts.imagePullSecret

        // -- labels
        if( opts.labels ) {
            def keys = opts.labels.keySet()
            if( 'app' in keys ) throw new IllegalArgumentException("Invalid pod label -- `app` is a reserved label")
            if( 'runName' in keys ) throw new IllegalArgumentException("Invalid pod label -- `runName` is a reserved label")
            labels.putAll( opts.labels )
        }

        // -- configMaps
        if( opts.getMountConfigMaps() )
            configMaps.addAll( opts.getMountConfigMaps() )

        // -- secrets
        if( opts.getMountSecrets() )
            secrets.addAll( opts.getMountSecrets() )

        // -- volume claims
        if( opts.getVolumeClaims() )
            volumeClaims.addAll( opts.getVolumeClaims() )

        // -- node selector
        if( opts.nodeSelector )
            nodeSelector = opts.nodeSelector

        // -- priority class name
        priorityClassName = opts.priorityClassName

        // -- security context
        if( opts.securityContext )
            securityContext = opts.securityContext

        return this
    }

    @PackageScope List<Map> createImagePullSecret() {
        def result = new ArrayList(1)
        def entry = new LinkedHashMap(1)
        entry.name = imagePullSecret
        result.add(entry)
        return result
    }

    Map build() {
        assert this.podName, 'Missing pod name'
        assert this.imageName, 'Missing pod container image'
        assert this.command, 'Missing pod command'

        final metadata = [
            name: podName,
            namespace: this.namespace ?: 'default'
        ]

        final container = [
            name: this.podName,
            image: this.imageName,
            command: this.command
        ]

        final spec = [
            containers: [ container ],
            restartPolicy: this.restartPolicy ?: 'Never',
        ]

        // affinity
        if( this.affinity )
            spec.affinity = this.affinity

        // annotations
        if( this.annotations )
            metadata.annotations = sanitize0(this.annotations, 'annotation')

        // automount service account token
        if( ! this.automountServiceAccountToken )
            spec.automountServiceAccountToken = false

        // environment variables
        final env = []
        for( PodEnv envVar : this.envVars ) {
            env.add(envVar.toSpec())
        }

        if ( env )
            container.env = env

        // image pull policy
        if( this.imagePullPolicy )
            container.imagePullPolicy = this.imagePullPolicy

        // image pull secret
        if( this.imagePullSecret )
            spec.imagePullSecrets = this.createImagePullSecret()

        // labels
        final labels = this.labels ?: [:]
        if( labels )
            metadata.labels = sanitize0(labels, 'label')

        // node selector
        if( this.nodeSelector )
            spec.nodeSelector = this.nodeSelector.toSpec()

        // priority class name
        if( this.priorityClassName )
            spec.priorityClassName = this.priorityClassName

        // resources
        final res = [:]
        if( this.cpus )
            res.cpu = this.cpus
        if( this.memory )
            res.memory = this.memory

        if( res ) {
            container.resources = [
                requests: res,
                limits: new HashMap<>(res)
            ]
        }

        // security context
        if( this.securityContext )
            spec.securityContext = this.securityContext.toSpec()

        // service account
        if( this.serviceAccount )
            spec.serviceAccountName = this.serviceAccount

        // work dir
        if( this.workDir )
            container.workingDir = this.workDir

        // add accelerator resource
        if( this.accelerator ) {
            container.resources = addAcceleratorResources(this.accelerator, container.resources as Map)
        }

        // add volumes
        final volumeMounts = []
        final volumes = []
        final namesMap = [:]

        // create a volume name for each unique volume claim
        for( String claimName : volumeClaims.collect { it.claimName }.unique() ) {
            final volumeName = nextVolumeName()
            namesMap[claimName] = volumeName
            volumes << [name: volumeName, persistentVolumeClaim: [claimName: claimName]]
        }

        // -- volume claims
        for( PodVolumeClaim entry : volumeClaims ) {
            // check if we already have a volume for the pvc
            final name = namesMap.get(entry.claimName)
            final claim = [name: name, mountPath: entry.mountPath ]
            if( entry.subPath )
                claim.subPath = entry.subPath
            if( entry.readOnly )
                claim.readOnly = entry.readOnly
            volumeMounts << claim
        }

        // -- configMap volumes
        for( PodMountConfig entry : configMaps ) {
            final name = nextVolumeName()
            configMapToSpec(name, entry, volumeMounts, volumes)
        }

        // -- hostPath volumes
        for( PodMountHostPath entry : hostPaths ) {
            final name = nextVolumeName()
            volumeMounts << [name: name, mountPath: entry.mountPath]
            volumes << [name: name, hostPath: [path: entry.hostPath]]
        }

        // -- secret volumes
        for( PodMountSecret entry : secrets ) {
            final name = nextVolumeName()
            secretToSpec(name, entry, volumeMounts, volumes)
        }

        if( volumeMounts )
            container.volumeMounts = volumeMounts
        if( volumes )
            spec.volumes = volumes

        // define top-level pod spec
        final pod = [
            apiVersion: 'v1',
            kind: 'Pod',
            metadata: metadata,
            spec: spec
        ]

        return pod
    }


    @PackageScope
    @CompileDynamic
    String getAcceleratorType(AcceleratorResource accelerator) {

        def type = accelerator.type ?: 'nvidia.com'

        if ( type.contains('/') )
            // Assume the user has fully specified the resource type.
            return type

        // Assume we're using GPU and update as necessary.
        if( !type.contains('.') ) type += '.com'
        type += '/gpu'

        return type
    }


    @PackageScope
    @CompileDynamic
    Map addAcceleratorResources(AcceleratorResource accelerator, Map res) {

        if( res == null )
            res = new LinkedHashMap(2)

        def type = getAcceleratorType(accelerator)

        if( accelerator.request ) {
            final req = res.requests ?: new LinkedHashMap<>(2)
            req.put(type, accelerator.request)
            res.requests = req
        }
        if( accelerator.limit ) {
            final lim = res.limits ?: new LinkedHashMap<>(2)
            lim.put(type, accelerator.limit)
            res.limits = lim
        }

        return res
    }

    @PackageScope
    @CompileDynamic
    static void secretToSpec(String volumeName, PodMountSecret entry, List volumeMounts, List volumes ) {
        assert entry

        final secret = [secretName: entry.secretName]
        if( entry.secretKey ) {
            secret.items = [ [key: entry.secretKey, path: entry.fileName ] ]
        }

        volumeMounts << [name: volumeName, mountPath: entry.mountPath]
        volumes << [name: volumeName, secret: secret ]
    }

    @PackageScope
    @CompileDynamic
    static void configMapToSpec(String volumeName, PodMountConfig entry, List<Map> volumeMounts, List<Map> volumes ) {
        assert entry

        final config = [name: entry.configName]
        if( entry.configKey ) {
            config.items = [ [key: entry.configKey, path: entry.fileName ] ]
        }

        volumeMounts << [name: volumeName, mountPath: entry.mountPath]
        volumes << [name: volumeName, configMap: config ]
    }

    protected Map sanitize0(Map map, String kind) {
        final result = new HashMap(map.size())
        for( Map.Entry entry : map )
            result.put(entry.key, sanitize0(entry.key, entry.value, kind))
        return result
    }

    /**
     * Valid label must be an empty string or consist of alphanumeric characters, '-', '_' or '.',
     * and must start and end with an alphanumeric character.
     *
     * @param value
     * @return
     */
    protected String sanitize0( key, value, String kind ) {
        def str = String.valueOf(value)
        if( str.length() > 63 ) {
            log.debug "K8s $kind exceeds allowed size: 63 -- offending name=$key value=$str"
            str = str.substring(0,63)
        }
        str = str.replaceAll(/[^a-zA-Z0-9\.\_\-]+/, '_')
        str = str.replaceAll(/^[^a-zA-Z0-9]+/, '')
        str = str.replaceAll(/[^a-zA-Z0-9]+$/, '')
        return str
    }

}
