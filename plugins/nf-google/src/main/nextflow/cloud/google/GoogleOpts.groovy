/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cloud.google

import com.google.auth.oauth2.GoogleCredentials
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.SysEnv
import nextflow.cloud.google.batch.client.BatchConfig
import nextflow.cloud.google.config.GoogleStorageOpts
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.exception.AbortOperationException
import nextflow.util.Duration
/**
 * Model Google config options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("google")
@Description("""
    The `google` scope allows you to configure the interactions with Google Cloud, including Google Cloud Batch and Google Cloud Storage.
""")
@Slf4j
@ToString(includeNames = true, includePackage = false)
@CompileStatic
class GoogleOpts implements ConfigScope {

    static final public String DEFAULT_LOCATION = 'us-central1'

    static Map<String,String> env = SysEnv.get()

    @ConfigOption
    @Description("""
        The Google Cloud project ID to use for pipeline execution.
    """)
    String project

    @ConfigOption
    @Description("""
        The Google Cloud location where jobs are executed (default: `us-central1`).
    """)
    final String location

    @ConfigOption
    @Description("""
        Use the given Google Cloud project ID as the billing project for storage access (default: `false`). Required when accessing data from [requester pays](https://cloud.google.com/storage/docs/requester-pays) buckets.
    """)
    final boolean enableRequesterPaysBuckets

    @ConfigOption
    @Description("""
        The HTTP connection timeout for Cloud Storage API requests (default: `'60s'`).
    """)
    final Duration httpConnectTimeout

    @ConfigOption
    @Description("""
        The HTTP read timeout for Cloud Storage API requests (default: `'60s'`).
    """)
    final Duration httpReadTimeout

    final BatchConfig batch

    final GoogleStorageOpts storage

    private File credsFile

    String getProjectId() { project }
    GoogleStorageOpts getStorageOpts() { storage }

    /* required by extension point -- do not remove */
    GoogleOpts() {}

    GoogleOpts(Map opts) {
        project = opts.project
        location = opts.location ?: DEFAULT_LOCATION
        enableRequesterPaysBuckets = opts.enableRequesterPaysBuckets as boolean
        httpConnectTimeout = opts.httpConnectTimeout ? opts.httpConnectTimeout as Duration : Duration.of('60s')
        httpReadTimeout = opts.httpReadTimeout ? opts.httpReadTimeout as Duration : Duration.of('60s')
        batch = new BatchConfig( opts.batch as Map ?: Collections.emptyMap() )
        storage = new GoogleStorageOpts( opts.storage as Map ?: Collections.emptyMap() )
    }

    @Memoized
    static GoogleOpts fromSession(Session session) {
        try {
            return fromSession0(session.config)
        }
        catch (Exception e) {
            if(session) session.abort()
            throw e
        }
    }

    protected static GoogleOpts fromSession0(Map config) {
        final result = new GoogleOpts( config.google as Map ?: Collections.emptyMap() )

        if( result.enableRequesterPaysBuckets && !result.projectId )
            throw new IllegalArgumentException("Config option 'google.enableRequesterPaysBuckets' cannot be honoured because the Google project Id has not been specified - Provide it by adding the option 'google.project' in the nextflow.config file")

        return result
    }

    static protected String getProjectIdFromCreds(String credsFilePath) {
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
        catch (Exception e) {
            throw new AbortOperationException("Invalid or corrupted Google credentials file: $credsFilePath", e)
        }
    }

    static GoogleOpts create(Session session) {
        // Typically the credentials picked up are the "Application Default Credentials"
        // as described at:
        //   https://github.com/googleapis/google-auth-library-java
        //
        // In that case, the project ID needs to be set in the nextflow config file.

        // If instead, the GOOGLE_APPLICATION_CREDENTIALS environment variable is set,
        // then the project ID will be picked up (along with the credentials) from the
        // JSON file that environment variable points to.

        final config = fromSession(session)

        def projectId
        def credsPath = env.get('GOOGLE_APPLICATION_CREDENTIALS')
        if( credsPath && (projectId = getProjectIdFromCreds(credsPath)) ) {
            config.credsFile = new File(credsPath)
            if( !config.project )
                config.project = projectId
            else if( config.project != projectId )
                throw new AbortOperationException("Project Id `$config.project` declared in the nextflow config file does not match the one expected by credentials file: $credsPath")
        }

        if( !config.project ) {
            throw new AbortOperationException("Missing Google project Id -- Specify it adding the setting `google.project='your-project-id'` in the nextflow.config file")
        }

        return config
    }

    @Memoized // make memoized to prevent multiple access to the creds file
    GoogleCredentials getCredentials() {
        return makeCreds(credsFile)
    }

    static protected GoogleCredentials makeCreds(File credsFile) {
        GoogleCredentials result
        if( credsFile ) {
            log.debug "Google auth via application credentials file: $credsFile"
            result = GoogleCredentials.fromStream(new FileInputStream(credsFile))
        }
        else {
            log.debug "Google auth via application DEFAULT"
            result = GoogleCredentials.getApplicationDefault()
        }
        return result.createScoped("https://www.googleapis.com/auth/cloud-platform")
    }
}
