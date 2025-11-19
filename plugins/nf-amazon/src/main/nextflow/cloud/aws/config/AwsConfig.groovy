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
import nextflow.config.spec.PlaceholderName
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

    final AwsS3ClientConfig client

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

    @PlaceholderName("<name>")
    final Map<String, AwsBucketConfig> buckets

    /* required by extension point -- do not remove */
    AwsConfig() {}

    AwsConfig(Map opts) {
        this.accessKey = opts.accessKey
        this.secretKey = opts.secretKey
        this.profile = getAwsProfile0(SysEnv.get(), opts)
        this.region = getAwsRegion(SysEnv.get(), opts)
        this.batch = new AwsBatchConfig((Map)opts.batch ?: Collections.emptyMap())
        this.client = new AwsS3ClientConfig((Map)opts.client ?: Collections.emptyMap())
        this.buckets = parseBuckets((Map<String,Map>)opts.buckets ?:Collections.<String,Map> emptyMap())
    }

    private static Map<String, AwsBucketConfig> parseBuckets(Map<String,Map> buckets){
        final result = new LinkedHashMap<String,AwsBucketConfig>()
        buckets.each { Map.Entry<String, Map> entry ->
            result[entry.key] = new AwsBucketConfig(entry.value)
        }
        return result
    }

    List<String> getCredentials() {
        return accessKey && secretKey
                ? List.of(accessKey, secretKey)
                : Collections.<String>emptyList()
    }

    AwsS3ClientConfig getS3Config() { client }

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
     *      1. bucket specific endpoint region
     *      2. global s3 client endpoint region
     *      3. bucket specific region
     *      5. config region
     *      6. US_EAST_1
     *
     *  @returns Resolved region.
     **/
    String resolveS3Region(String bucketName) {
        def bucketRegion = null
        def bucketEpRegion = null
        if( bucketName && buckets && buckets.containsKey(bucketName) ){
            bucketRegion = buckets[bucketName].region
            bucketEpRegion =  buckets[bucketName].getEndpointRegion()
        }
        return bucketEpRegion ?: client.getEndpointRegion() ?: bucketRegion ?: this.region ?: Region.US_EAST_1.id()
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

    Map getS3LegacyProperties(String bucketName) {
        final result = new LinkedHashMap(20)

        // -- global client config options
        def config = client.getAwsClientConfig()

        // -- overwrite with bucket specific options
        if( bucketName && buckets && buckets.containsKey(bucketName) ){
            config = config + buckets[bucketName].toLegacyConfig()
        }

        config = checkDefaultErrorRetry(config, SysEnv.get())
        if( config ) {
            result.putAll(config)
        }

        log.debug "AWS S3 config properties: ${dumpAwsConfig(result)}"
        return result
    }

    AwsBucketConfig getBucketConfig(String bucketName){
        // Get global bucket
        def config = client.toBucketConfigMap()

        // overwrite with bucket specific options
        if( bucketName && buckets && buckets.containsKey(bucketName) ){
            config = config + buckets[bucketName].toBucketConfigMap()
        }

        return new AwsBucketConfig(config)
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

    String generateUploadCliArgs(String bucketName){
        final config = getBucketConfig(bucketName)
        final region = config.region ?: region
        final cliArgs = []
        cliArgs.add(config.storageClass ? "--storage-class ${config.storageClass}" : "--storage-class STANDARD")
        if(!batch.s5cmdPath && region)
            cliArgs.add("--region ${region}")
        if( config.anonymous )
            cliArgs.add("--no-sign-request")
        if( config.storageEncryption )
            cliArgs.add("--sse ${config.storageEncryption}")
        if( config.storageKmsKeyId )
            cliArgs.add("--sse-kms-key-id ${config.storageKmsKeyId}")
        if( config.s3Acl )
            cliArgs.add("--acl ${config.s3Acl}")
        if( config.requesterPays )
            cliArgs.add("--request-payer requester")
        if( config.endpoint)
            cliArgs.add("--endpoint-url ${config.endpoint}")
        return cliArgs.join(" ")
    }

    String generateDownloadCliArgs(String bucketName){
        final config = getBucketConfig(bucketName)
        final region = config.region ?: region
        final cliArgs = []
        if(!batch.s5cmdPath && region)
            cliArgs.add("--region ${region}")
        if( config.anonymous )
            cliArgs.add("--no-sign-request")
        if( config.endpoint)
            cliArgs.add("--endpoint-url ${config.endpoint}")
        return cliArgs.join(" ")
    }
}
