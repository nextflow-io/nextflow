/*
 * Copyright 2013-2024, Seqera Labs
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


import java.util.function.Consumer

import groovy.util.logging.Slf4j
import nextflow.util.Duration
import nextflow.util.IniFile
import nextflow.util.MemoryUnit
import nextflow.util.TestOnly
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

    @TestOnly
    static void reset() {
        session = null
        config = null
        hooks.clear()
    }
}
