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

package nextflow.cloud

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CloudDriverFactoryTest extends Specification {

    def 'should load at least on driver' () {

        when:
        def found = CloudDriverFactory.loadDrivers()
        then:
        found.size()>0
        found.fake == FakeCloudDriver.class

    }

    def 'should return a set of names'() {
        expect:
        CloudDriverFactory.getDriverNames()  == ['fake','aws'] as Set
    }

    def 'should return the driver instance' () {
        expect:
        CloudDriverFactory.getDriver('fake') instanceof FakeCloudDriver
    }
}
