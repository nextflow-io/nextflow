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

import nextflow.Const
import org.gridgain.grid.GridConfiguration
import org.gridgain.grid.cache.GridCacheConfiguration
import org.gridgain.grid.logger.slf4j.GridSlf4jLogger

/**
 * Creates a the GridGain configuration object common to master and worker nodes
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GgConfig {

    static final SESSIONS_CACHE = 'allSessions'

    static final GRID_NAME = Const.APP_NAME

    static GridConfiguration create( String role ) {
        assert role in ['master','worker'], "Parameter 'role' can be either 'master' or 'worker'"

        // configure
        System.setProperty('GRIDGAIN_DEBUG_ENABLED','true')
        System.setProperty('GRIDGAIN_UPDATE_NOTIFIER','false')

        def cfg = new GridConfiguration()
        cfg.setGridName( GRID_NAME )
        cfg.setUserAttributes( ROLE: role )
        cfg.setGridLogger( new GridSlf4jLogger() )
        //cfg.setDeploymentSpi( new GgDeploymentSpi() )
        //cfg.setPeerClassLoadingEnabled(true)
        //cfg.setDiscoverySpi( new GridTcpDiscoverySpi() )

        // configure session cache
        def sessionCfg = new GridCacheConfiguration()
        sessionCfg.setName(SESSIONS_CACHE)
        sessionCfg.setStartSize(128)
        cfg.setCacheConfiguration( sessionCfg )

        return cfg
    }

}
