/*
 * Copyright 2018, WuxiNextcode
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
package nextflow.cloud.google.pipelines

import spock.lang.Specification

class GooglePipelinesConfigurationTest extends Specification {

    def 'should construct correctly'() {
        given:
        def name = "testName"
        def testZone = ["testZone1","testZone2"]
        def testRegion = ["region1,region2"]
        def preemp = true
        def remoteBinDir = null

        when:
        def config1 = new GooglePipelinesConfiguration(name,testZone,testRegion,remoteBinDir,preemp)
        def config2 = new GooglePipelinesConfiguration(name,testZone,testRegion)

        then:
        with(config1) {
            project == name
            zone == testZone
            region == testRegion
            preemptible == preemp
        }

        with(config2) {
            project == name
            zone == testZone
            region == testRegion
            !preemptible
        }
    }
}
