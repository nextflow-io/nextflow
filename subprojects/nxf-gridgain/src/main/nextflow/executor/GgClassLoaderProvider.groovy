/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import groovy.transform.Memoized
import nextflow.util.RemoteSession
import org.gridgain.grid.Grid
import org.gridgain.grid.cache.GridCache
import org.gridgain.grid.logger.GridLogger
import org.gridgain.grid.resources.GridInstanceResource
import org.gridgain.grid.resources.GridLoggerResource

/**
 * Creates a class-loader based on the client classpath.
 *
 * @see GgClosureTask
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GgClassLoaderProvider {

    @GridLoggerResource
    private GridLogger log;


    @GridInstanceResource
    private Grid grid;


    @Memoized
    ClassLoader getClassLoaderFor( UUID sessionId ) {
        assert sessionId
        GridCache<UUID, RemoteSession> allSessions = grid.cache( GgGridFactory.SESSIONS_CACHE )

        if( !allSessions )
            throw new IllegalStateException('Missing session cache object')

        def session = allSessions.get(sessionId)
        if( !session )
            throw new IllegalStateException("Missing session object for id: $sessionId")

        def loader = new GroovyClassLoader()
        session.classpath.each { File file ->
            log.debug "Adding to classpath: $file"
            loader.addClasspath(file.absolutePath)
        }

        return loader
    }

    @Memoized
    RemoteSession getSessionFor( UUID sessionId ) {
        assert sessionId
        GridCache<UUID, RemoteSession> allSessions = grid.cache( GgGridFactory.SESSIONS_CACHE )

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
