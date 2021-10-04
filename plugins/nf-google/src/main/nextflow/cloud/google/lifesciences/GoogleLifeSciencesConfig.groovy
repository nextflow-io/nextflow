/*
 * Copyright 2019, Google Inc
 * Copyright 2018, WuxiNextcode
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

package nextflow.cloud.google.lifesciences

import java.nio.file.Path

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.cloud.CloudTransferOptions
import nextflow.exception.AbortOperationException
import nextflow.util.Duration
import nextflow.util.MemoryUnit

/**
 * Helper class wrapping configuration required for Google Pipelines.
 *
 * @author Ólafur Haukur Flygenring <olafurh@wuxinextcode.com>
 * @author  Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString(includePackage = false, includeNames = true)
@CompileStatic
class GoogleLifeSciencesConfig implements CloudTransferOptions {

    public final static String DEFAULT_COPY_IMAGE = 'google/cloud-sdk:slim'

    public final static String DEFAULT_SSH_IMAGE = 'gcr.io/cloud-genomics-pipelines/tools'

    public final static String DEFAULT_ENTRY_POINT = '/bin/bash'

    public final static int DEF_PARALLEL_THREAD_COUNT =1
    public final static int DEF_DOWNLOAD_MAX_COMPONENTS =8

    String project
    List<String> zones
    List<String> regions
    boolean preemptible
    Path remoteBinDir
    String location
    boolean disableBinDir
    MemoryUnit bootDiskSize
    String cpuPlatform
    boolean sshDaemon
    String sshImage
    Integer debugMode
    String copyImage
    boolean usePrivateAddress
    boolean enableRequesterPaysBuckets
    String network
    String subnetwork
    String serviceAccountEmail
    int parallelThreadCount
    int downloadMaxComponents
    boolean keepAliveOnFailure

    int maxParallelTransfers = MAX_TRANSFER
    int maxTransferAttempts = MAX_TRANSFER_ATTEMPTS
    Duration delayBetweenAttempts = DEFAULT_DELAY_BETWEEN_ATTEMPTS

    @Deprecated
    GoogleLifeSciencesConfig(String project, List<String> zone, List<String> region, Path remoteBinDir = null, boolean preemptible = false) {
        this.project = project
        this.zones = zone
        this.regions = region
        this.remoteBinDir = remoteBinDir
        this.preemptible = preemptible
        // infer location
        this.location = region ? region.get(0) : null
        if( !location )
            throw new IllegalArgumentException("Missing Google cloud location")
    }

    GoogleLifeSciencesConfig() {}

    @Memoized
    static GoogleLifeSciencesConfig fromSession(Session session) {
        try {
            fromSession0(session.config)
        }
        catch (Exception e) {
            if(session) session.abort()
            throw e
        }
    }

    protected static GoogleLifeSciencesConfig fromSession0(Map config) {
        def project = config.navigate("google.project") as String

        //check if we have one of the mutual exclusive zone or region specified
        if(!config.navigate("google.zone") && !config.navigate("google.region")){
            throw new AbortOperationException("Missing configuration value 'google.zone' or 'google.region'")
        }

        //check if we have one of the mutual exclusive zone or region specified
        if(config.navigate("google.zone") && config.navigate("google.region")){
            throw new AbortOperationException("You can't specify both 'google.zone' and 'google.region' configuration parameters -- Please remove one of them from your configuration")
        }

        def path = config.navigate('env.PATH')
        if( path ) {
            log.warn "Environment PATH defined in config file is ignored by Google Pipeline executor"
        }

        /*
         * upload local binaries
         */
        final boolean disableBinDir = config.navigate('google.lifeSciences.disableRemoteBinDir',false)
        final preemptible = config.navigate("google.lifeSciences.preemptible", false) as boolean
        final bootDiskSize = config.navigate('google.lifeSciences.bootDiskSize') as MemoryUnit
        final cpuPlatform = config.navigate('google.lifeSciences.cpuPlatform') as String
        final sshDaemon = config.navigate('google.lifeSciences.sshDaemon', false) as boolean
        final sshImage = config.navigate('google.lifeSciences.sshImage', DEFAULT_SSH_IMAGE) as String
        final copyImage = config.navigate('google.lifeSciences.copyImage', DEFAULT_COPY_IMAGE) as String
        final debugMode = config.navigate('google.lifeSciences.debug', System.getenv('NXF_DEBUG'))
        final privateAddr  = config.navigate('google.lifeSciences.usePrivateAddress') as boolean
        final requesterPays = config.navigate('google.enableRequesterPaysBuckets') as boolean

        final network = config.navigate('google.lifeSciences.network') as String
        final subnetwork = config.navigate('google.lifeSciences.subnetwork') as String
        final serviceAccountEmail = config.navigate('google.lifeSciences.serviceAccountEmail') as String
        final keepAlive = config.navigate('google.lifeSciences.keepAliveOnFailure', false) as boolean

        def zones = (config.navigate("google.zone") as String)?.split(",")?.toList() ?: Collections.<String>emptyList()
        def regions = (config.navigate("google.region") as String)?.split(",")?.toList() ?: Collections.<String>emptyList()
        def location = config.navigate("google.location") as String ?: fallbackToRegionOrZone(regions,zones)

        final maxParallelTransfers = config.navigate('google.storage.maxParallelTransfers', MAX_TRANSFER) as int
        final maxTransferAttempts = config.navigate('google.storage.maxTransferAttempts', MAX_TRANSFER_ATTEMPTS) as int
        final delayBetweenAttempts = config.navigate('google.storage.delayBetweenAttempts', DEFAULT_DELAY_BETWEEN_ATTEMPTS) as Duration
        final parallelThreadCount = config.navigate('google.storage.parallelThreadCount', DEF_PARALLEL_THREAD_COUNT) as int
        final downloadMaxComponents = config.navigate('google.storage.downloadMaxComponents', DEF_DOWNLOAD_MAX_COMPONENTS) as int


        new GoogleLifeSciencesConfig(
                project: project,
                regions: regions,
                zones: zones,
                location: location,
                preemptible: preemptible,
                disableBinDir: disableBinDir,
                bootDiskSize: bootDiskSize,
                cpuPlatform: cpuPlatform,
                debugMode: debugMode0(debugMode),
                copyImage: copyImage,
                sshDaemon: sshDaemon,
                sshImage: sshImage,
                usePrivateAddress: privateAddr,
                enableRequesterPaysBuckets: requesterPays,
                maxParallelTransfers: maxParallelTransfers,
                maxTransferAttempts: maxTransferAttempts,
                delayBetweenAttempts: delayBetweenAttempts,
                network: network,
                subnetwork: subnetwork,
                serviceAccountEmail: serviceAccountEmail,
                parallelThreadCount: parallelThreadCount,
                downloadMaxComponents: downloadMaxComponents,
                keepAliveOnFailure: keepAlive
            )
    }

    static private Integer debugMode0(value) {
        if( value instanceof Boolean )
            return value ? 1 : 0
        else if( value ) {
            return value.toString().toInteger()
        }
        return null
    }

    static String fallbackToRegionOrZone(List<String> regions, List<String> zones) {
        if( regions ) {
            return bestLocationForRegion(regions[0])
        }
        if( zones ) {
            def norm = zones
                    .collect { int p = zones[0].lastIndexOf('-'); p!=-1 ? it.substring(0,p) : it }
                    .unique()
            return bestLocationForRegion(norm[0])
        }
        throw new AbortOperationException("Missing Google region or zone information")
    }

    static String bestLocationForRegion(String region) {
        region.startsWith('europe-') ? 'europe-west2' : 'us-central1'
    }

    static String getProjectIdFromCreds(String credsFilePath) {
        if( !credsFilePath )
            throw new AbortOperationException('Missing Google credentials -- make sure your environment defines the GOOGLE_APPLICATION_CREDENTIALS environment variable')

        final file = new File(credsFilePath)
        try {
            final creds = (Map)new JsonSlurper().parse(file)
            if( creds.project_id )
                return creds.project_id
            else
                throw new AbortOperationException("Missing `project_id` in Google credentials file: $credsFilePath")
        }
        catch(FileNotFoundException e) {
            throw new AbortOperationException("Missing Google credentials file: $credsFilePath")
        }
    }
}
