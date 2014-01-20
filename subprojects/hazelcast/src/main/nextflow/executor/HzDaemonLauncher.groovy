/*
 * Copyright (c) 2012, the authors.
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

import com.hazelcast.config.Config
import com.hazelcast.core.Hazelcast
import com.hazelcast.core.HazelcastInstance
import groovy.util.logging.Slf4j
import nextflow.daemon.DaemonLauncher

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class HzDaemonLauncher implements DaemonLauncher {

    HazelcastInstance hazelcast

    HzDaemonLauncher() { }

    @Override
    void launch(Map properties) {
        log.debug "Launching Hazelcast instance"
        def cfg = new Config();
        cfg.setProperty('hazelcast.logging.type', 'slf4j')
        cfg.setProperty('hazelcast.system.log.enabled','true')
        hazelcast = Hazelcast.newHazelcastInstance(cfg)
    }
}
