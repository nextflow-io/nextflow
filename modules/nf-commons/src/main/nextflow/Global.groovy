/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow

import java.nio.file.Path
import java.nio.file.Paths
import java.util.function.Consumer

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.util.Duration
import nextflow.util.IniFile
import nextflow.util.MemoryUnit
import org.apache.commons.lang.StringUtils
import org.apache.commons.lang.exception.ExceptionUtils
/**
 * Hold global variables
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class Global {

    /**
     * The pipeline session instance
     */
    static private ISession session

    /**
     * Allow to load session in a lazy manner
     */
    static private Closure<ISession> loader

    /**
     * The main configuration object
     */
    static Map config

    /**
     * @return The object instance representing the current session
     */
    static ISession getSession() {
        if( session != null )
            return session
        if( loader != null )
            session = loader.call()
        return session
    }

    /**
     * Set the application session object
     *
     * @param value An object instance representing the current session
     */
    static void setSession( ISession value ) {
        session = value
    }

    /**
     * Set a session lazy loader
     *
     * @param loader
     */
    static void setLazySession( Closure<ISession> loader ) {
        Global.loader = loader
    }

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
    static List<String> getAwsCredentials0( Map env, Map config, List<Path> files = []) {

        String a
        String b

        if( config && config.aws instanceof Map ) {
            a = ((Map)config.aws).accessKey
            b = ((Map)config.aws).secretKey

            if( a && b ) {
                log.debug "Using AWS credentials defined in nextflow config file"
                return [a, b]
            }

        }

        // as define by amazon doc
        // http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html
        if( env && (a=env.AWS_ACCESS_KEY_ID) && (b=env.AWS_SECRET_ACCESS_KEY) )  {
            log.debug "Using AWS credentials defined by environment variables AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY"
            return [a, b]
        }

        if( env && (a=env.AWS_ACCESS_KEY) && (b=env.AWS_SECRET_KEY) ) {
            log.debug "Using AWS credentials defined by environment variables AWS_ACCESS_KEY and AWS_SECRET_KEY"
            return [a, b]
        }

        for( Path it : files ) {
            final conf = new IniFile(it)
            final profile = getAwsProfile0(env, config)
            final section = conf.section(profile)
            if( (a=section.aws_access_key_id) && (b=section.aws_secret_access_key) ) {
                final token = section.aws_session_token
                if( token ) {
                    log.debug "Using AWS temporary session credentials defined in `$profile` section in file: ${conf.file}"
                    return [a,b,token]
                }
                else {
                    log.debug "Using AWS credential defined in `$profile` section in file: ${conf.file}"
                    return [a,b]
                }
            }
        }

        return null
    }

    static protected String getAwsProfile0(Map env, Map<String,Object> config) {

        final profile = config?.navigate('aws.profile')
        if( profile )
            return profile

        if( env?.containsKey('AWS_PROFILE'))
            return env.get('AWS_PROFILE')

        if( env?.containsKey('AWS_DEFAULT_PROFILE'))
            return env.get('AWS_DEFAULT_PROFILE')

        return 'default'
    }

    static List<String> getAwsCredentials(Map env, Map config) {

        def home = Paths.get(System.properties.get('user.home') as String)
        def files = [ home.resolve('.aws/credentials'), home.resolve('.aws/config') ]
        getAwsCredentials0(env, config, files)

    }

    static String getAwsRegion(Map env=null, Map config=null) {
        if( env==null ) env = SysEnv.get()
        if( config==null ) config = this.config

        def home = Paths.get(System.properties.get('user.home') as String)
        def file = home.resolve('.aws/config')

        return getAwsRegion0(env, config, file)
    }

    static protected String getAwsRegion0(Map env, Map config, Path awsFile) {
        // check nxf config file
        if( config && config.aws instanceof Map ) {
            def region = ((Map)config.aws).region
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

    static List<String> getAwsCredentials(Map env) {
        getAwsCredentials(env, config)
    }

    static List<String> getAwsCredentials() {
        getAwsCredentials(SysEnv.get(), config)
    }

    static Map<String,?> getAwsClientConfig() {
        if( config?.aws?.client instanceof Map ) {
            return normalizeAwsClientConfig(config.aws.client as Map)
        }

        return null
    }

    static String getAwsS3Endpoint() {
        getAwsS3Endpoint0(SysEnv.get(), config ?: Collections.emptyMap())
    }

    static protected String getAwsS3Endpoint0(Map<String,String> env, Map<String,Object> config) {
        config.navigate('aws.client.endpoint', env.get('AWS_S3_ENDPOINT'))
    }

    /**
     * Convert configuration keys from camel-case notation (nextflow) to underscore
     * separated notation expected by the AWS client
     *
     * @return A map object containing the AWS client configuration properties
     */
    static protected Map normalizeAwsClientConfig(Map<String,?> client) {

        normalizeMemUnit(client, 'uploadChunkSize');
        normalizeDuration(client, 'uploadRetrySleep');


        def result = [:]
        client.each { String name, value ->
            def newKey = name.isCamelCase() ? StringUtils.splitByCharacterTypeCamelCase(name).join('_').toLowerCase() : name
            result.put(newKey,value?.toString())
        }
        return result
    }

    static void normalizeMemUnit(Map<String,?> client, String key) {
        if( client.get(key) instanceof String ) {
            client.put(key, MemoryUnit.of((String)client.get(key)))
        }
        if( client.get(key) instanceof MemoryUnit ) {
            client.put(key, ((MemoryUnit)client.get(key)).toBytes())
        }
    }

    static void normalizeDuration(Map<String,?> client, String key)  {
        if( client.get(key) instanceof String ) {
            client.put(key, Duration.of((String)client.get(key)))
        }
        if( client.get(key) instanceof Duration ) {
            client.put(key, ((Duration)client.get(key)).toMillis())
        }
    }

    /**
     * Run the specified closure on application termination
     *
     * @param callback A closure to be executed on application shutdown
     */
    static void onCleanup(Consumer<ISession> callback) {
        if( callback==null ) {
            log.warn "Cleanup consumer cannot be null\n${ExceptionUtils.getStackTrace(new Exception())}"
            return 
        }
        hooks.add(callback)
    }

    static final private List<Consumer<ISession>> hooks = []

    static synchronized cleanUp() {
        for( Consumer<ISession> c : hooks ) {
            try {
                c.accept(session)
            }
            catch( Exception e ) {
                log.debug("Error during on cleanup", e )
            }
        }
    }

}
