/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 * Hold global variables
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class Global {

    /**
     * The pipeline session instance
     */
    static private ISession session

    /**
     * The main configuration object
     */
    static Map config

    /**
     * @return The object instance representing the current session
     */
    static ISession getSession() {
        session
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
    static List<String> getAwsCredentials( Map env, Map config ) {

        String a
        String b

        if( config && config.aws instanceof Map ) {
            a = ((Map)config.aws).accessKey
            b = ((Map)config.aws).secretKey

            if( a && b )
                return [a, b]
        }

        if( env && (a=env.AWS_ACCESS_KEY) && (b=env.AWS_SECRET_KEY) ) {
            return [a, b]
        }

        // as define by amazon doc
        // http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html
        if( env && (a=env.AWS_ACCESS_KEY_ID) && (b=env.AWS_SECRET_ACCESS_KEY) )  {
            return [a, b]
        }

        return null
    }

    static List<String> getAwsCredentials(Map env) {
        getAwsCredentials(env, config)
    }

    static List<String> getAwsCredentials() {
        getAwsCredentials(System.getenv(), config)
    }

    /**
     * Run the specified closure on application termination
     *
     * @param callback A closure to be executed on application shutdown
     */
    static void onShutdown(Closure callback) {
        hooks.add(callback)
    }

    static final List<Closure> hooks = []

    static private synchronized cleanUp() {
        for( Closure c : hooks ) {
            try {
                c.call()
            }
            catch( Exception e ) {
                log.debug("Error during on cleanup", e )
            }
        }
    }

    /*
     * Global shutdown hook
     */
    static {
        Runtime.getRuntime().addShutdownHook {
            cleanUp()
        }
    }


}
