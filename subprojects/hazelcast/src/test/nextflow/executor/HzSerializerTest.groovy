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

import java.nio.file.Paths

import com.hazelcast.config.Config
import com.hazelcast.core.Hazelcast
import com.hazelcast.core.HazelcastInstance
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HzSerializerTest extends Specification {

    def testGlobalSerializer() {

        setup:
        Config cfg = new Config()
        HzSerializerConfig.registerAll( cfg.getSerializationConfig() )
        HazelcastInstance hz = Hazelcast.newHazelcastInstance(cfg)
        Map map = hz.getMap('map')
        final path = Paths.get('/some/file')

        when:
        map.put( 'path', path )
        then:
        map.get('path') == Paths.get('/some/file')


        when:
        def uuid = UUID.randomUUID()
        def obj = [x:1, path: path, uuid: uuid ]
        map.put('item', obj)
        then:
        map.get('item') == obj

        when:
        def closure = { 'Hello' }
        closure.delegate = [x:1, path: Paths.get('/hola')]
        def cmd = new HzCmdCall(UUID.randomUUID(),1, path, closure )
        map.put( 'cmd', cmd )
        def result = new HzCmdResult(cmd, 'Bravo', new IllegalArgumentException() )
        map.put( 'result', result )
        then:
        map.get('cmd') == cmd
        map.get('result').sender == result.sender
        map.get('result').context == result.context
        map.get('result').taskId == result.taskId

    }


}
