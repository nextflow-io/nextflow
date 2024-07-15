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
import nextflow.exception.ProcessUnrecoverableException
import nextflow.util.MemoryUnit
/**
 * Model Google Batch config settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BatchConfig {

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
    private BatchRetryConfig retryConfig
    private List<Integer> autoRetryExitCodes

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
    BatchRetryConfig getRetryConfig() { retryConfig }
    List<Integer> getAutoRetryExitCodes() { autoRetryExitCodes }

    static BatchConfig create(Session session) {
        final result = new BatchConfig()
        result.googleOpts = GoogleOpts.create(session)
        result.credentials = result.googleOpts.credentials
        result.allowedLocations = session.config.navigate('google.batch.allowedLocations', List.of()) as List<String>
        result.bootDiskSize = session.config.navigate('google.batch.bootDiskSize') as MemoryUnit
        result.cpuPlatform = session.config.navigate('google.batch.cpuPlatform')
        result.maxSpotAttempts = session.config.navigate('google.batch.maxSpotAttempts',5) as int
        result.installGpuDrivers = session.config.navigate('google.batch.installGpuDrivers',false)
        result.preemptible = session.config.navigate('google.batch.preemptible',false)
        result.spot = session.config.navigate('google.batch.spot',false)
        result.usePrivateAddress = session.config.navigate('google.batch.usePrivateAddress',false)
        result.network = session.config.navigate('google.batch.network')
        result.subnetwork = session.config.navigate('google.batch.subnetwork')
        result.serviceAccountEmail = session.config.navigate('google.batch.serviceAccountEmail')
        result.retryConfig = new BatchRetryConfig( session.config.navigate('google.batch.retryPolicy') as Map ?: Map.of() )
        result.autoRetryExitCodes = parseAutoRetryExitCodes0(session.config.navigate('google.batch.autoRetryExitCodes','50001'))
        return result
    }

    @Override
    String toString(){
        return "BatchConfig[googleOpts=$googleOpts"
    }

    static private String _50001 = '50001'

    static private List<Integer> DEFAULT_RETRY_LIST = List.of(_50001.toInteger())

    static protected List<Integer> parseAutoRetryExitCodes0(value) {
        if(!value)
            return List.of()
        if( value instanceof List ) {
            // it's expected to be a list of integer
            return value as List<Integer>
        }
        if( value instanceof CharSequence ) {
            final v0 = value.toString()
            if( v0==_50001 )
                return DEFAULT_RETRY_LIST
            else
                return v0.tokenize(',').toList().collect(it->parseCode(it))
        }
        return null
    }

    static private Integer parseCode(String v) {
        try {
            return v.toInteger()
        }
        catch (NumberFormatException e) {
            throw new ProcessUnrecoverableException("Invalid exit code value: $v -- check the setting `google.batch.autoRetryExitCodes`")
        }
    }

}
