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

    static enum MetaType { LABEL, ANNOTATION }

    static enum SegmentType {
        PREFIX (253),
        NAME (63),
        VALUE (63)

        private final int maxSize;
        SegmentType(int maxSize) {
            this.maxSize = maxSize;
        }
    }

    static @PackageScope AtomicInteger VOLUMES = new AtomicInteger()

    String podName

    String imageName

    String imagePullPolicy

    String imagePullSecret

    List<String> command = []

    List<String> args = new ArrayList<>()

    Map<String,String> labels = [:]

    Map<String,String> annotations = [:]

    String namespace

    String restart

    List<PodEnv> envVars = []

    String workDir

    Integer cpus

    String memory

    String disk

    String serviceAccount

    boolean automountServiceAccountToken = true

    AcceleratorResource accelerator

    Collection<PodMountConfig> configMaps = []

    Collection<PodMountCsiEphemeral> csiEphemerals = []

    Collection<PodMountEmptyDir> emptyDirs = []

    Collection<PodMountSecret> secrets = []

    Collection<PodHostMount> hostMounts = []

    Collection<PodVolumeClaim> volumeClaims = []

    PodSecurityContext securityContext

    PodNodeSelector nodeSelector

    Map affinity

    String priorityClassName

    List<Map> tolerations = []

    boolean privileged

    int activeDeadlineSeconds

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
        if( cmd==null ) return this
        assert cmd instanceof List || cmd instanceof CharSequence, "Missing or invalid K8s command parameter: $cmd"
        this.command = cmd instanceof List ? cmd : ['/bin/bash','-c', cmd.toString()]
        return this
    }

    PodSpecBuilder withArgs( args ) {
        if( args==null ) return this
        assert args instanceof List || args instanceof CharSequence, "Missing or invalid K8s args parameter: $args"
        this.args = args instanceof List ? args : ['/bin/bash','-c', args.toString()]
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

    PodSpecBuilder withDisk(String disk) {
        this.disk = disk
        return this
    }

    PodSpecBuilder withDisk(MemoryUnit disk)  {
        this.disk = "${disk.mega}Mi".toString()
        return this
    }

    PodSpecBuilder withAccelerator(AcceleratorResource acc) {
        this.accelerator = acc
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

    PodSpecBuilder withAnnotation( String name, String value ) {
        this.annotations.put(name, value)
        return this
    }

    PodSpecBuilder withAnnotations(Map annotations) {
        this.annotations.putAll(annotations)
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

    PodSpecBuilder withCsiEphemerals( Collection<PodMountCsiEphemeral> csiEphemerals ) {
        this.csiEphemerals.addAll(csiEphemerals)
        return this
    }

    PodSpecBuilder withCsiEphemeral( PodMountCsiEphemeral csiEphemeral ) {
        this.csiEphemerals.add(csiEphemeral)
        return this
    }

    PodSpecBuilder withEmptyDirs( Collection<PodMountEmptyDir> emptyDirs ) {
        this.emptyDirs.addAll(emptyDirs)
        return this
    }

    PodSpecBuilder withEmptyDir( PodMountEmptyDir emptyDir ) {
        this.emptyDirs.add(emptyDir)
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

    PodSpecBuilder withPrivileged(boolean value) {
        this.privileged = value
        return this
    }

    PodSpecBuilder withActiveDeadline(int seconds) {
        this.activeDeadlineSeconds = seconds
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
        // -- configMaps
        if( opts.getMountConfigMaps() )
            configMaps.addAll( opts.getMountConfigMaps() )
        // -- csi ephemeral volumes
        if( opts.getMountCsiEphemerals() )
            csiEphemerals.addAll( opts.getMountCsiEphemerals() )
        // -- emptyDirs
        if( opts.getMountEmptyDirs() )
            emptyDirs.addAll( opts.getMountEmptyDirs() )
        // -- secrets
        if( opts.getMountSecrets() )
            secrets.addAll( opts.getMountSecrets() )
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
        // - annotations
        if( opts.annotations ) {
            annotations.putAll( opts.annotations )
        }
        // -- security context
        if( opts.securityContext )
            securityContext = opts.securityContext
        // -- node selector
        if( opts.nodeSelector )
            nodeSelector = opts.nodeSelector
        // -- affinity
        if( opts.affinity )
            affinity = opts.affinity
        // -- automount service account token
        automountServiceAccountToken = opts.automountServiceAccountToken
        // -- priority class name
        priorityClassName = opts.priorityClassName
        // -- tolerations
        if( opts.tolerations )
            tolerations.addAll(opts.tolerations)
        // -- privileged
        privileged = opts.privileged

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
        assert this.command || this.args, 'Missing K8s command parameter'

        final restart = this.restart ?: 'Never'

        final metadata = new LinkedHashMap<String,Object>()
        metadata.name = podName
        metadata.namespace = namespace ?: 'default'

        final labels = this.labels ?: [:]
        final env = []
        for( PodEnv entry : this.envVars ) {
            env.add(entry.toSpec())
        }

        final container = [ name: this.podName, image: this.imageName ]
        if( this.command )
            container.command = this.command
        if( this.args )
            container.args = args

        if( this.workDir )
            container.put('workingDir', workDir)

        if( imagePullPolicy )
            container.imagePullPolicy = imagePullPolicy

        if( privileged ) {
            // note: privileged flag needs to be defined in the *container* securityContext
            // not the 'spec' securityContext (see below)
            container.securityContext = [ privileged: true ]
        }

        final spec = [
                restartPolicy: restart,
                containers: [ container ],
        ]

        if( nodeSelector )
            spec.nodeSelector = nodeSelector.toSpec()

        if( affinity )
            spec.affinity = affinity

        if( this.serviceAccount )
            spec.serviceAccountName = this.serviceAccount

        if( ! this.automountServiceAccountToken )
            spec.automountServiceAccountToken = false

        if( securityContext )
            spec.securityContext = securityContext.toSpec()

        if( imagePullSecret )
            spec.imagePullSecrets = createPullSecret()

        if( priorityClassName )
            spec.priorityClassName = priorityClassName

        // tolerations
        if( this.tolerations )
            spec.tolerations = this.tolerations

        // add labels
        if( labels )
            metadata.labels = sanitize(labels, MetaType.LABEL)

        if( annotations )
            metadata.annotations = sanitize(annotations, MetaType.ANNOTATION)

        // time directive
        if ( activeDeadlineSeconds > 0)
            spec.activeDeadlineSeconds = activeDeadlineSeconds

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
        if( this.cpus ) {
            container.resources = addCpuResources(this.cpus, container.resources as Map)
        }

        if( this.memory ) {
            container.resources = addMemoryResources(this.memory, container.resources as Map)
        }

        if( this.accelerator ) {
            container.resources = addAcceleratorResources(this.accelerator, container.resources as Map)
        }

        if( this.disk ) {
            container.resources = addDiskResources(this.disk, container.resources as Map)
        }

        // add storage definitions ie. volumes and mounts
        final List<Map> mounts = []
        final List<Map> volumes = []
        final namesMap = [:]

        // creates a volume name for each unique claim name
        for( String claimName : volumeClaims.collect { it.claimName }.unique() ) {
            final volName = nextVolName()
            namesMap[claimName] = volName
            volumes << [name: volName, persistentVolumeClaim: [claimName: claimName]]
        }

        // -- persistent volume claims
        for( PodVolumeClaim entry : volumeClaims ) {
            //check if we already have a volume for the pvc
            final name = namesMap.get(entry.claimName)
            final claim = [name: name, mountPath: entry.mountPath ]
            if( entry.subPath )
                claim.subPath = entry.subPath
            if( entry.readOnly )
                claim.readOnly = entry.readOnly
            mounts << claim
        }

        // -- configMap volumes
        for( PodMountConfig entry : configMaps ) {
            final name = nextVolName()
            configMapToSpec(name, entry, mounts, volumes)
        }

        // -- csi ephemeral volumes
        for( PodMountCsiEphemeral entry : csiEphemerals ) {
            final name = nextVolName()
            mounts << [name: name, mountPath: entry.mountPath, readOnly: entry.csi.readOnly ?: false]
            volumes << [name: name, csi: entry.csi]
        }

        // -- emptyDir volumes
        for( PodMountEmptyDir entry : emptyDirs ) {
            final name = nextVolName()
            mounts << [name: name, mountPath: entry.mountPath]
            volumes << [name: name, emptyDir: entry.emptyDir]
        }

        // -- secret volumes
        for( PodMountSecret entry : secrets ) {
            final name = nextVolName()
            secretToSpec(name, entry, mounts, volumes)
        }

        // -- host path volumes
        for( PodHostMount entry : hostMounts ) {
            final name = nextVolName()
            mounts << [name: name, mountPath: entry.mountPath]
            volumes << [name: name, hostPath: [path: entry.hostPath]]
        }


        if( volumes )
            spec.volumes = volumes
        if( mounts )
            container.volumeMounts = mounts

        return pod
    }

    Map buildAsJob() {
        final pod = build()

        // job metadata
        final metadata = new LinkedHashMap<String,Object>()
        metadata.name = this.podName    //  just use the podName for simplicity, it may be renamed to just `name` or `resourceName` in the future
        metadata.namespace = this.namespace ?: 'default'

        // job spec
        final spec = new LinkedHashMap<String,Object>()
        spec.backoffLimit = 0
        spec.template = [spec: pod.spec]

        if( labels )
            metadata.labels = sanitize(labels, MetaType.LABEL)

        if( annotations )
            metadata.annotations = sanitize(annotations, MetaType.ANNOTATION)

        final result = [
                apiVersion: 'batch/v1',
                kind: 'Job',
                metadata: metadata,
                spec: spec ]

        return result

    }

    @PackageScope
    Map addCpuResources(Integer cpus, Map res) {
        if( res == null )
            res = [:]

        final req = res.requests as Map ?: new LinkedHashMap<>(10)
        req.cpu = cpus
        res.requests = req

        return res
    }

    @PackageScope
    Map addMemoryResources(String memory, Map res) {
        if( res == null )
            res = new LinkedHashMap(10)

        final req = res.requests as Map ?: new LinkedHashMap(10)
        req.memory = memory
        res.requests = req

        final lim = res.limits as Map ?: new LinkedHashMap(10)
        lim.memory = memory
        res.limits = lim

        return res
    }

    @PackageScope
    Map addDiskResources(String diskRequest, Map res) {
        if( res == null )
            res = new LinkedHashMap(10)

        final req = res.requests as Map ?: new LinkedHashMap(10)
        req.'ephemeral-storage' = diskRequest
        res.requests = req

        final lim = res.limits as Map ?: new LinkedHashMap(10)
        lim.'ephemeral-storage' = diskRequest
        res.limits = lim

        return res
    }

    @PackageScope
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
    Map addAcceleratorResources(AcceleratorResource accelerator, Map res) {

        if( res == null )
            res = new LinkedHashMap(2)

        def type = getAcceleratorType(accelerator)

        if( accelerator.request ) {
            final req = res.requests as Map ?: new LinkedHashMap<>(2)
            req.put(type, accelerator.request)
            res.requests = req
        }
        if( accelerator.limit ) {
            final lim = res.limits as Map ?: new LinkedHashMap<>(2)
            lim.put(type, accelerator.limit)
            res.limits = lim
        }

        return res
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

    protected Map sanitize(Map map, MetaType kind) {
        final result = new HashMap(map.size())
        for( Map.Entry entry : map ) {
            final key = sanitizeKey(entry.key as String, kind)
            final value = (kind == MetaType.LABEL)
                ? sanitizeValue(entry.value, kind, SegmentType.VALUE)
                : entry.value

            result.put(key, value)
        }
        return result
    }

    protected String sanitizeKey(String value, MetaType kind) {
        final parts = value.tokenize('/')
        
        if (parts.size() == 2) {
            return "${sanitizeValue(parts[0], kind, SegmentType.PREFIX)}/${sanitizeValue(parts[1], kind, SegmentType.NAME)}"
        }
        if( parts.size() == 1 ) {
            return sanitizeValue(parts[0], kind, SegmentType.NAME)
        }
        else {
            throw new IllegalArgumentException("Invalid key in pod ${kind.toString().toLowerCase()} -- Key can only contain exactly one '/' character")
        }
    }


    /**
     * Sanitize a string value to contain only alphanumeric characters, '-', '_' or '.',
     * and to start and end with an alphanumeric character.
     */
    protected String sanitizeValue(value, MetaType kind, SegmentType segment) {
        def str = String.valueOf(value)
        if( str.length() > segment.maxSize ) {
            log.debug "K8s $kind $segment exceeds allowed size: $segment.maxSize -- offending str=$str"
            str = str.substring(0,segment.maxSize)
        }
        str = str.replaceAll(/[^a-zA-Z0-9\.\_\-]+/, '_')
        str = str.replaceAll(/^[^a-zA-Z0-9]+/, '')
        str = str.replaceAll(/[^a-zA-Z0-9]+$/, '')
        return str
    }

}
