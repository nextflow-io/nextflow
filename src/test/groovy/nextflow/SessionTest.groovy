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

package nextflow

import nextflow.util.Duration
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SessionTest extends Specification {


    def testBaseDirAndBinDir() {

        setup:
        def base = File.createTempDir()
        def bin = new File(base,'bin'); bin.mkdir()

        when:
        def session = new Session()
        then:
        session.baseDir == null
        session.binDir == null

        when:
        session = new Session()
        session.baseDir = new File('some/folder')
        then:
        session.baseDir == new File('some/folder')
        session.binDir == null

        when:
        session = new Session()
        session.baseDir = base
        then:
        session.baseDir == base
        session.binDir.exists()

        cleanup:
        base.deleteDir()

    }


    def testGetQueueSize() {

        when:
        def session = [:] as Session
        session.config = [ executor:['$sge':[queueSize: 123] ] ]
        then:
        session.getQueueSize('sge', 1) == 123
        session.getQueueSize('xxx', 1) == 1
        session.getQueueSize(null, 1) == 1

        when:
        def session2 = [:] as Session
        session2.config = [ executor:[ queueSize: 321, '$sge':[queueSize:789] ] ]
        then:
        session2.getQueueSize('sge', 2) == 789
        session2.getQueueSize('xxx', 2) == 321
        session2.getQueueSize(null, 2) == 321


        when:
        def session3 = [:] as Session
        session3.config = [ executor: 'sge' ]
        then:
        session3.getQueueSize('sge', 1) == 1
        session3.getQueueSize('xxx', 2) == 2
        session3.getQueueSize(null, 3) == 3


    }

    def testGetPollInterval() {

        when:
        def session1 = [:] as Session
        session1.config = [ executor:['$sge':[pollInterval: 345] ] ]
        then:
        session1.getPollInterval('sge').toMillis() == 345
        session1.getPollInterval('xxx').toMillis() == 1_000
        session1.getPollInterval(null).toMillis() == 1_000
        session1.getPollInterval(null, 2_000 as Duration).toMillis() == 2_000

        when:
        def session2 = [:] as Session
        session2.config = [ executor:[ pollInterval: 321, '$sge':[pollInterval:789] ] ]
        then:
        session2.getPollInterval('sge').toMillis() == 789
        session2.getPollInterval('xxx').toMillis() == 321
        session2.getPollInterval(null).toMillis() == 321

        when:
        def session3 = [:] as Session
        session3.config = [ executor: 'lsf' ]
        then:
        session3.getPollInterval('sge', 33 as Duration ).toMillis() == 33
        session3.getPollInterval('xxx', 44 as Duration ).toMillis() == 44
        session3.getPollInterval(null, 55 as Duration).toMillis() == 55

    }

    def testGetExitReadTimeout() {

        setup:
        def session1 = [:] as Session
        session1.config = [ executor:['$sge':[exitReadTimeout: '5s'] ] ]

        expect:
        session1.getExitReadTimeout('sge') == '5sec' as Duration
        session1.getExitReadTimeout('lsf', '3sec' as Duration) == '3sec' as Duration

    }

    def testGetQueueStatInterval() {

        setup:
        def session1 = [:] as Session
        session1.config = [ executor:['$sge':[queueStatInterval: '4sec'] ] ]

        expect:
        session1.getQueueStatInterval('sge') == '4sec' as Duration
        session1.getQueueStatInterval('lsf', '1sec' as Duration) == '1sec' as Duration

    }

    def testMonitorDumpInterval() {

        setup:
        def session1 = [:] as Session
        session1.config = [ executor:['$sge':[dumpInterval: '6sec'] ] ]

        expect:
        session1.getMonitorDumpInterval('sge') == '6sec' as Duration
        session1.getMonitorDumpInterval('lsf', '2sec' as Duration) == '2sec' as Duration

    }

    def testGetExecConfigProp() {

        when:
        def session = [:] as Session
        session.config = [ executor: [x:123, y:222, '$hazelcast': [y:333] ] ]
        then:
        session.getExecConfigProp( 'hazelcast', 'x', null ) == 123
        session.getExecConfigProp( 'hazelcast', 'y', null ) == 333
        session.getExecConfigProp( 'local', 'y', null ) == 222
        session.getExecConfigProp( 'local', 'y', 'beta') == 222
        session.getExecConfigProp( 'hazelcast', 'z', null ) ==  null
        session.getExecConfigProp( 'hazelcast', 'z', 'alpha') == 'alpha'
        session.getExecConfigProp( 'hazelcast', 'z', 'alpha', [NXF_EXECUTOR_Z:'hola']) == 'hola'
        session.getExecConfigProp( 'hazelcast', 'p.q.z', null, [NXF_EXECUTOR_P_Q_Z:'hello']) == 'hello'
    }

}
