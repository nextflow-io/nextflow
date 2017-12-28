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

package nextflow.util

import spock.lang.IgnoreIf
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BarrierTest extends Specification {

    def 'test is terminated' () {

        when:
        def barrier = new Barrier()
        then:
        barrier.isCompleted()

        when:
        barrier.register(1)
        then:
        !barrier.isCompleted()

        when:
        barrier.register(1)
        barrier.register(2)
        barrier.arrive(1)
        then:
        !barrier.isCompleted()

        when:
        barrier.arrive(2)
        then:
        barrier.isCompleted()

    }

    @IgnoreIf({ javaVersion == 1.7 })
    def 'test await termination' () {

        given:
        def barrier = new Barrier()
        def begin = System.currentTimeMillis()

        barrier.register(1)
        barrier.register(2)
        Thread.start { sleep 50; barrier.arrive(2) }
        Thread.start { sleep 100; barrier.arrive(1) }

        when:
        barrier.awaitCompletion()
        def elapsed = System.currentTimeMillis() - begin

        then:
        elapsed>=100 && elapsed<500
    }
}
