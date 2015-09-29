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

package nextflow.executor
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.util.RemoteSession
import org.apache.ignite.Ignite
import org.apache.ignite.IgniteCache
import org.apache.ignite.IgniteLogger
import org.apache.ignite.resources.IgniteInstanceResource
import org.apache.ignite.resources.LoggerResource

/**
 * Creates a class-loader based on the client classpath.
 *
 * @see IgClosureTask
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class IgClassLoaderProvider {

    @LoggerResource
    private IgniteLogger log;

    @IgniteInstanceResource
    private Ignite grid;

    @Memoized
    ClassLoader getClassLoaderFor( UUID sessionId ) {
        assert sessionId
        def allSessions = (IgniteCache<UUID, RemoteSession>)grid.cache( IgGridFactory.SESSIONS_CACHE )

        if( !allSessions )
            throw new IllegalStateException('Missing session cache object')

        def session = allSessions.get(sessionId)
        if( !session )
            throw new IllegalStateException("Missing session object for id: $sessionId")

        def loader = new GroovyClassLoader()
        session.classpath.each { Path file ->
            log.debug "Adding to classpath: $file"
            loader.addClasspath(file.toAbsolutePath().toString())
        }

        return loader
    }

    @Memoized
    RemoteSession getSessionFor( UUID sessionId ) {
        assert sessionId
        IgniteCache<UUID, RemoteSession> allSessions = grid.cache( IgGridFactory.SESSIONS_CACHE )

        if( !allSessions )
            throw new IllegalStateException('Missing session cache object')

        def session = allSessions.get(sessionId)
        if( !session )
            throw new IllegalStateException("Missing session object for id: $sessionId")

        return session
    }

    @Memoized
    Map getConfigFor( UUID sessionId ) {
        getSessionFor(sessionId).getConfig()
    }

}
