/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
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

    GoogleOpts getGoogleOpts() { return googleOpts }
    GoogleCredentials getCredentials() { return credentials }
    boolean getDisableBinDir() { disableBinDir }

    static BatchConfig create(Session session) {
        final result = new BatchConfig()
        result.googleOpts = GoogleOpts.create(session)
        result.credentials = makeCreds(result.googleOpts.credsFile)
        result.disableBinDir = session.config.navigate('google.batch.disableRemoteBinDir',false)
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
