/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
