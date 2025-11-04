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

import software.amazon.awssdk.regions.Region
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.SysEnv
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.IniFile
/**
 * Model AWS cloud configuration settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("aws")
@Description("""
    The `aws` scope controls the interactions with AWS, including AWS Batch and S3.
""")
@Slf4j
@CompileStatic
class AwsConfig implements ConfigScope {

    final AwsBatchConfig batch

    final AwsS3Config client

    @ConfigOption
    @Description("""
        AWS region (e.g. `us-east-1`).
    """)
    final String region

    @ConfigOption
    @Description("""
        AWS account access key.
    """)
    final String accessKey

    @ConfigOption
    @Description("""
        AWS account secret key.
    """)
    final String secretKey

    @ConfigOption
    @Description("""
        AWS profile from `~/.aws/credentials`.
    """)
    final String profile

    /* required by extension point -- do not remove */
    AwsConfig() {}

    AwsConfig(Map opts) {
        this.accessKey = opts.accessKey
        this.secretKey = opts.secretKey
        this.profile = getAwsProfile0(SysEnv.get(), opts)
        this.region = getAwsRegion(SysEnv.get(), opts)
        this.batch = new AwsBatchConfig((Map)opts.batch ?: Collections.emptyMap())
        this.client = new AwsS3Config((Map)opts.client ?: Collections.emptyMap())
    }

    List<String> getCredentials() {
        return accessKey && secretKey
                ? List.of(accessKey, secretKey)
                : Collections.<String>emptyList()
    }

    AwsS3Config getS3Config() { client }

    AwsBatchConfig getBatchConfig() { batch }

    @Deprecated
    String getS3GlobalRegion() {
        return !region || !s3Config.endpoint || s3Config.endpoint.contains(".amazonaws.com")
            ? Region.US_EAST_1.id()         // always use US_EAST_1 as global region for AWS endpoints
            : region                        // for custom endpoint use the config provided region
    }

    /**
     *  Resolves the region used for S3 evaluating the region resolved from config and a possible region defined in the endpoint.
     *  Fallback to the global region US_EAST_1 when no region is found.
     *
     *  Preference:
     *      1. endpoint region
     *      2. config region
     *      3. US_EAST_1
     *
     *  @returns Resolved region.
     **/
    String resolveS3Region() {
        final epRegion = client.getEndpointRegion()
        return epRegion ?: this.region ?: Region.US_EAST_1.id()
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
        if( env && env.AWS_REGION )  {
            return env.AWS_REGION.toString()
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
        def config = client.getAwsClientConfig()
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
