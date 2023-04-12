/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.cloud.aws.config

import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.SysEnv
import nextflow.util.IniFile
/**
 * Model AWS cloud configuration settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsConfig {

    private AwsBatchConfig batchConfig

    private AwsS3Config s3config

    private String region

    private String accessKey

    private String secretKey

    private String profile
    
    private AwsS3Legacy s3Legacy

    AwsConfig(Map config) {
        this.accessKey = config.accessKey
        this.secretKey = config.secretKey
        this.profile = getAwsProfile0(SysEnv.get(), config)
        this.region = getAwsRegion(SysEnv.get(), config)
        this.batchConfig = new AwsBatchConfig((Map)config.batch ?: Collections.emptyMap())
        this.s3config = new AwsS3Config((Map)config.client ?: Collections.emptyMap())
        this.s3Legacy = new AwsS3Legacy((Map)config.client ?: Collections.emptyMap())
    }

    String getAccessKey() { accessKey }

    String getSecretKey() { secretKey }

    List<String> getCredentials() {
        return accessKey && secretKey
                ? List.of(accessKey, secretKey)
                : Collections.<String>emptyList()
    }

    String getProfile() { profile }

    String getRegion() { region }

    AwsS3Config getS3Config() { s3config }

    AwsBatchConfig getBatchConfig() { batchConfig }

    Map<String,?> getS3LegacyClientConfig() {
        return s3Legacy.getAwsClientConfig()
    }

    static protected String getAwsProfile0(Map env, Map<String,Object> config) {

        final profile = config?.profile as String
        if( profile )
            return profile

        if( env?.containsKey('AWS_PROFILE'))
            return env.get('AWS_PROFILE')

        if( env?.containsKey('AWS_DEFAULT_PROFILE'))
            return env.get('AWS_DEFAULT_PROFILE')

        return null
    }


    static protected String getAwsRegion(Map env, Map config) {

        def home = Paths.get(System.properties.get('user.home') as String)
        def file = home.resolve('.aws/config')

        return getAwsRegion0(env, config, file)
    }

    static protected String getAwsRegion0(Map env, Map config, Path awsFile) {
        // check nxf config file
        if( config instanceof Map ) {
            def region = config.region
            if( region )
                return region.toString()
        }

        if( env && env.AWS_DEFAULT_REGION )  {
            return env.AWS_DEFAULT_REGION.toString()
        }

        if( !awsFile.exists() ) {
            return null
        }

        final profile = getAwsProfile0(env, config) ?: 'default'
        final ini = new IniFile(awsFile)
        return ini.section(profile).region
    }

    Map getS3LegacyProperties() {
        final result = new LinkedHashMap(20)

        // -- remaining client config options
        def config = getS3LegacyClientConfig()
        config = checkDefaultErrorRetry(config, SysEnv.get())
        if( config ) {
            result.putAll(config)
        }

        log.debug "AWS S3 config properties: ${dumpAwsConfig(result)}"
        return result
    }

    static protected Map checkDefaultErrorRetry(Map result, Map env) {
        if( result == null )
            result = new HashMap(10)

        if( result.max_error_retry==null ) {
            result.max_error_retry = env?.AWS_MAX_ATTEMPTS
        }
        // fallback to default
        if( result.max_error_retry==null ) {
            result.max_error_retry = '5'
        }
        // make sure that's a string value as it's expected by the client
        else {
            result.max_error_retry = result.max_error_retry.toString()
        }

        return result
    }

    static private String dumpAwsConfig( Map<String,String> config ) {
        def result = new HashMap(config)
        if( config.access_key && config.access_key.size()>6 )
            result.access_key = "${config.access_key.substring(0,6)}.."

        if( config.secret_key && config.secret_key.size()>6 )
            result.secret_key = "${config.secret_key.substring(0,6)}.."

        if( config.session_token && config.session_token.size()>6 )
            result.session_token = "${config.session_token.substring(0,6)}.."

        return result.toString()
    }

    static private AwsConfig getConfig0(Map config) {
        if( config==null ) {
            log.warn("Missing nextflow session config object")
            return new AwsConfig(Collections.emptyMap())
        }
        new AwsConfig( (Map)config.aws ?: Collections.emptyMap()  )
    }

    static AwsConfig config() {
        getConfig0(Global.config)
    }
}
