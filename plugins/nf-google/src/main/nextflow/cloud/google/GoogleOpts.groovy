/*
 * Copyright 2020-2021, Seqera Labs
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
 *
 */

package nextflow.cloud.google

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.ToString
import nextflow.Session
import nextflow.exception.AbortOperationException
/**
 * Model Google config options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@CompileStatic
class GoogleOpts {

    static Map<String,String> env = System.getenv()

    private String projectId
    private List<String> zones
    private List<String> regions
    private String location
    private File credsFile
    private boolean enableRequesterPaysBuckets

    String getProjectId() { projectId }
    File getCredsFile() { credsFile }
    String getLocation() { location }
    boolean getEnableRequesterPaysBuckets() { enableRequesterPaysBuckets }

    @Memoized
    static GoogleOpts fromSession(Session session) {
        try {
            fromSession0(session.config)
        }
        catch (Exception e) {
            if(session) session.abort()
            throw e
        }
    }

    protected static GoogleOpts fromSession0(Map config) {
        final result = new GoogleOpts()
        result.projectId = config.navigate("google.project") as String

        //check if we have one of the mutual exclusive zone or region specified
        if(!config.navigate("google.zone") && !config.navigate("google.region")){
            throw new AbortOperationException("Missing configuration value 'google.zone' or 'google.region'")
        }

        //check if we have one of the mutual exclusive zone or region specified
        if(config.navigate("google.zone") && config.navigate("google.region")){
            throw new AbortOperationException("You can't specify both 'google.zone' and 'google.region' configuration parameters -- Please remove one of them from your configuration")
        }

        result.zones = (config.navigate("google.zone") as String)?.split(",")?.toList() ?: Collections.<String>emptyList()
        result.regions = (config.navigate("google.region") as String)?.split(",")?.toList() ?: Collections.<String>emptyList()
        result.location = config.navigate("google.location") as String ?: fallbackToRegionOrZone(result.regions, result.zones)
        result.enableRequesterPaysBuckets = config.navigate('google.enableRequesterPaysBuckets') as boolean

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
            if( !config.projectId )
                config.projectId = projectId
            else if( config.projectId != projectId )
                throw new AbortOperationException("Project Id `$config.projectId` declared in the nextflow config file does not match the one expected by credentials file: $credsPath")
        }

        if( !config.projectId ) {
            throw new AbortOperationException("Missing Google project Id -- Specify it adding the setting `google.project='your-project-id'` in the nextflow.config file")
        }

        return config
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
}
