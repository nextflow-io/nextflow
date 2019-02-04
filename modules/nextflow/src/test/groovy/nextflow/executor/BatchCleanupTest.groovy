/*
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
