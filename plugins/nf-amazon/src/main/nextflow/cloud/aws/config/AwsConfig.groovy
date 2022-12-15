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
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
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

    private AwsBatchConfig batchOpts

    private AwsS3Config s3Opts

    private String region

    private String accessKey

    private String secretKey

    private String profile
    
    private String assumeRoleArn

    AwsConfig(Map config) {
        final creds = getAwsCredentials(SysEnv.get(), config)
        this.accessKey = creds?[0]
        this.secretKey = creds?[1]
        this.profile = getAwsProfile0(SysEnv.get(), config)
        this.region = getAwsRegion(SysEnv.get(), config)
        this.assumeRoleArn = config.assumeRoleArn as String
        this.batchOpts = new AwsBatchConfig( (Map)config.batch ?: Collections.emptyMap() )
        this.s3Opts = new AwsS3Config( (Map)config.s3 ?: Collections.emptyMap() )
    }

    String getAccessKey() { accessKey }

    String getSecretKey() { secretKey }

    String getProfile() { profile }

    String getRegion() { region }

    String getAssumeRoleArn() { assumeRoleArn }

    /**
     * Retrieve the AWS credentials from the given context. It look for AWS credential in the following order
     * 1) Nextflow config {@code aws.accessKey} and {@code aws.secretKey} pair
     * 2) System env {@code AWS_ACCESS_KEY} and {@code AWS_SECRET_KEY} pair
     * 3) System env {@code AWS_ACCESS_KEY_ID} and {@code AWS_SECRET_ACCESS_KEY} pair
     *
     *
     * @param env The system environment map
     * @param config The nextflow config object map
     * @return A pair where the first element is the access key and the second the secret key or
     *      {@code null} if the credentials are missing
     */
    @PackageScope
    static List<String> getAwsCredentials0(Map env, Map config, List<Path> files=List.of()) {

        String a
        String b

        if( config instanceof Map ) {
            a = config.accessKey
            b = config.secretKey

            if( a && b ) {
                log.debug "Using AWS credentials defined in nextflow config file"
                return List.of(a, b)
            }

        }

        // as define by amazon doc
        // http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html
        if( env && (a=env.AWS_ACCESS_KEY_ID) && (b=env.AWS_SECRET_ACCESS_KEY) )  {
            log.debug "Using AWS credentials defined by environment variables AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY"
            return List.of(a,b)
        }

        if( env && (a=env.AWS_ACCESS_KEY) && (b=env.AWS_SECRET_KEY) ) {
            log.debug "Using AWS credentials defined by environment variables AWS_ACCESS_KEY and AWS_SECRET_KEY"
            return List.of(a, b)
        }

        for( Path it : files ) {
            final conf = new IniFile(it)
            final profile = getAwsProfile0(env, config)
            final section = conf.section(profile)
            if( (a=section.aws_access_key_id) && (b=section.aws_secret_access_key) ) {
                log.debug "Using AWS credential defined in `$profile` section in file: ${conf.file}"
                return List.of(a,b)
            }
        }

        return null
    }

    static protected String getAwsProfile0(Map env, Map<String,Object> config) {

        final profile = config?.profile as String
        if( profile )
            return profile

        if( env?.containsKey('AWS_PROFILE'))
            return env.get('AWS_PROFILE')

        if( env?.containsKey('AWS_DEFAULT_PROFILE'))
            return env.get('AWS_DEFAULT_PROFILE')

        return 'default'
    }

    static List<String> getAwsCredentials(Map env, Map config) {

        final home = Paths.get(System.properties.get('user.home') as String)
        final files = [ home.resolve('.aws/credentials'), home.resolve('.aws/config') ]
        return getAwsCredentials0(env, config, files)

    }

    static String getAwsRegion(Map env, Map config) {

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

        def profile = getAwsProfile0(env, config)
        def ini = new IniFile(awsFile)
        return ini.section(profile).region
    }

}
