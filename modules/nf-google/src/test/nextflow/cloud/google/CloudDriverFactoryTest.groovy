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

package nextflow.cloud.google

import spock.lang.Requires
import spock.lang.Specification

import nextflow.cloud.CloudDriverFactory

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CloudDriverFactoryTest extends Specification {


    def 'should return a set of names'() {
        expect:
        CloudDriverFactory.getDriverNames().contains('google')
    }

    @Requires({ System.getenv('GOOGLE_APPLICATION_CREDENTIALS') })
    def 'should return the driver instance' () {
        given:
        def cfg = [project:'rare-lattice-222412', zone:'europe-west1-b']
        expect:
        CloudDriverFactory.getDriver('google', cfg) instanceof GoogleCloudDriver
    }
}
