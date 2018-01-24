/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BatchCleanupTest extends Specification {

    def 'should collect jobs' () {

        given:
        def batch = new BatchCleanup()
        def lsf = Stub(LsfExecutor)
        lsf.getName() >> 'lsf'
        def sge = Stub(SgeExecutor)
        sge.getName() >> 'sge'

        when:
        batch.collect(lsf, 100)
        batch.collect(lsf, 101)
        batch.collect(sge, 311)

        then:
        batch.aggregate.size() == 2
        batch.aggregate.get('lsf').jobIds == [100,101]
        batch.aggregate.get('sge').jobIds == [311]

    }


    def 'should kill jobs in groups' () {

        given:
        def batch = new BatchCleanup()
        def lsf = Mock(AbstractGridExecutor)
        lsf.getName() >> 'lsf'
        and:
        batch.size = 5

        when:
        5.times { batch.collect(lsf, 100+it)  }
        and:
        batch.kill()
        then:
        1 * lsf.killTask( [100,101,102,103,104] )


        when:
        batch.aggregate.clear()
        and:
        10.times { batch.collect(lsf, 100+it)  }
        and:
        batch.collect(lsf, 201)
        batch.collect(lsf, 202)
        batch.collect(lsf, 203)
        and:
        batch.kill()

        then:
        1 * lsf.killTask( [100,101,102,103,104] )
        then:
        1 * lsf.killTask( [105,106,107,108,109] )
        then:
        1 * lsf.killTask( [201,202,203] )


    }

}
