/*
 * Copyright 2022, Google Inc.
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

package nextflow.cloud.google.batch.client


import com.google.auth.oauth2.GoogleCredentials
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.cloud.google.GoogleOpts
import nextflow.util.ConfigHelper
import nextflow.util.MemoryUnit
/**
 * Model Google Batch config settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BatchConfig {

    static private final Set<String> VALID_OPTIONS = [
        'allowedLocations',
        'bootDiskSize',
        'cpuPlatform',
        'network',
        'preemptible',
        'serviceAccountEmail',
        'spot',
        'subnetwork',
        'usePrivateAddress',
    ]

    private GoogleOpts googleOpts
    private GoogleCredentials credentials
    private List<String> allowedLocations
    private MemoryUnit bootDiskSize
    private String cpuPlatform
    private int maxSpotAttempts
    private boolean installGpuDrivers
    private boolean preemptible
    private boolean spot
    private boolean usePrivateAddress
    private String network
    private String subnetwork
    private String serviceAccountEmail

    GoogleOpts getGoogleOpts() { return googleOpts }
    GoogleCredentials getCredentials() { return credentials }
    List<String> getAllowedLocations() { allowedLocations }
    MemoryUnit getBootDiskSize() { bootDiskSize }
    String getCpuPlatform() { cpuPlatform }
    int getMaxSpotAttempts() { maxSpotAttempts }
    boolean getInstallGpuDrivers() { installGpuDrivers }
    boolean getPreemptible() { preemptible }
    boolean getSpot() { spot }
    boolean getUsePrivateAddress() { usePrivateAddress }
    String getNetwork() { network }
    String getSubnetwork() { subnetwork }
    String getServiceAccountEmail() { serviceAccountEmail }

    static BatchConfig create(Session session) {
        final config = session.config.navigate('google.batch', [:]) as Map
        ConfigHelper.checkInvalidConfigOptions('google.batch', config, VALID_OPTIONS)

        final result = new BatchConfig()
        result.googleOpts = GoogleOpts.create(session)
        result.credentials = result.googleOpts.credentials
        result.allowedLocations = (config.allowedLocations ?: List.of()) as List<String>
        result.bootDiskSize = config.bootDiskSize as MemoryUnit
        result.cpuPlatform = config.cpuPlatform
        result.maxSpotAttempts = (config.maxSpotAttempts ?: 5) as int
        result.installGpuDrivers = config.installGpuDrivers ?: false
        result.preemptible = config.preemptible ?: false
        result.spot = config.spot ?: false
        result.usePrivateAddress = config.usePrivateAddress ?: false
        result.network = config.network
        result.subnetwork = config.subnetwork
        result.serviceAccountEmail = config.serviceAccountEmail
        return result
    }

    @Override
    String toString(){
        return "BatchConfig[googleOpts=$googleOpts"
    }

}
