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
/**
 * Model Google Batch config settings
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BatchConfig {

    private GoogleOpts googleOpts
    private GoogleCredentials credentials
    private boolean disableBinDir
    private boolean spot
    private boolean preemptible
    private boolean usePrivateAddress
    private String network
    private String subnetwork
    private String serviceAccountEmail

    GoogleOpts getGoogleOpts() { return googleOpts }
    GoogleCredentials getCredentials() { return credentials }
    boolean getDisableBinDir() { disableBinDir }
    boolean getPreemptible() { preemptible }
    boolean getSpot() { spot }
    boolean getUsePrivateAddress() { usePrivateAddress }
    String getNetwork() { network }
    String getSubnetwork() { subnetwork }
    String getServiceAccountEmail() { serviceAccountEmail }

    static BatchConfig create(Session session) {
        final result = new BatchConfig()
        result.googleOpts = GoogleOpts.create(session)
        result.credentials = makeCreds(result.googleOpts.credsFile)
        result.disableBinDir = session.config.navigate('google.batch.disableRemoteBinDir',false)
        result.spot = session.config.navigate('google.batch.spot',false)
        result.preemptible = session.config.navigate('google.batch.preemptible',false)
        result.usePrivateAddress = session.config.navigate('google.batch.usePrivateAddress',false)
        result.network = session.config.navigate('google.batch.network')
        result.subnetwork = session.config.navigate('google.batch.subnetwork')
        result.serviceAccountEmail = session.config.navigate('google.batch.serviceAccountEmail')
        return result
    }

    static protected GoogleCredentials makeCreds(File credsFile) {
        GoogleCredentials result
        if( credsFile ) {
            log.debug "Google auth via application credentials file: $credsFile"
            result = GoogleCredentials .fromStream(new FileInputStream(credsFile))
        }
        else {
            log.debug "Google auth via application DEFAULT"
            result = GoogleCredentials.getApplicationDefault()
        }
        return result.createScoped("https://www.googleapis.com/auth/cloud-platform")
    }

    @Override
    String toString(){
        return "BatchConfig[googleOpts=$googleOpts; disableBinDir=$disableBinDir]"
    }

}
