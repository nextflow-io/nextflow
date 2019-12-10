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
package nextflow.cloud.google.lifesciences

import nextflow.exception.AbortOperationException
import spock.lang.Specification

class GoogleLifeSciencesConfigTest extends Specification {

    def 'should construct correctly'() {
        given:
        def name = "testName"
        def testZone = ["testZone1","testZone2"]
        def testRegion = ["region1,region2"]
        def preemp = true
        def remoteBinDir = null

        when:
        def config1 = new GoogleLifeSciencesConfig(project: name, zones: testZone, regions:testRegion,remoteBinDir: remoteBinDir,preemptible: preemp)
        def config2 = new GoogleLifeSciencesConfig(project: name, zones: testZone, regions:testRegion)

        then:
        with(config1) {
            project == name
            zones == testZone
            regions == testRegion
            preemptible == preemp
        }

        with(config2) {
            project == name
            zones == testZone
            regions == testRegion
            !preemptible
        }
    }

    def 'should config location ' () {
        given:
        def map = [google:
                           [
                                   project: 'foo',
                                   zone: 'us-east4-a,us-east4-b',
                                   location: 'eu-west1',
                           ]]

        when:
        def config = GoogleLifeSciencesConfig.fromSession0(map)
        then:
        config.project == 'foo'
        config.regions == []
        config.zones == ['us-east4-a','us-east4-b']
        config.location == 'eu-west1'
        config.preemptible
        !config.disableBinDir
    }

    def 'should config flags' () {
        when:
        def map = [google:
                           [
                                   project: 'bar',
                                   zone: 'eu-west1-a',
                                   lifeSciences: [preemptible: false, disableRemoteBinDir: false]
                           ]]
        def config = GoogleLifeSciencesConfig.fromSession0(map)
        then:
        config.project == 'bar'
        config.location == 'eu-west1'
        !config.preemptible
        !config.disableBinDir


        when:
        map = [google:
                           [
                                   project: 'bar',
                                   zone: 'eu-west1-a',
                                   lifeSciences: [preemptible: true, disableRemoteBinDir: true]
                           ]]
        config = GoogleLifeSciencesConfig.fromSession0(map)
        then:
        config.project == 'bar'
        config.location == 'eu-west1'
        config.preemptible
        config.disableBinDir
    }

    def 'should config location from region' () {
        given:
        def map = [google:
                           [ region: 'us-central', project: 'foo'] ]

        when:
        def config = GoogleLifeSciencesConfig.fromSession0(map)
        then:
        config.project == 'foo'
        config.regions == ['us-central']
        config.zones == []
        config.location == 'us-central'
    }

    def 'should config location from zone' () {
        given:
        def map = [google:
                           [ zone: 'us-east4-a,us-east4-c', project: 'foo'] ]

        when:
        def config = GoogleLifeSciencesConfig.fromSession0(map)
        then:
        config.project == 'foo'
        config.regions == []
        config.zones == ['us-east4-a','us-east4-c']
        config.location == 'us-east4'
    }

    def 'should report missing project' () {
        when:
        GoogleLifeSciencesConfig.fromSession0([:])
        then:
        def err =thrown(AbortOperationException)
        and:
        err.message.startsWith('Missing Google project Id')
    }

    def 'should report missing region' () {
        when:
        GoogleLifeSciencesConfig.fromSession0([google:[project:'foo']])
        then:
        def err =thrown(AbortOperationException)
        and:
        err.message.startsWith('Missing configuration value \'google.zone\' or \'google.region\'')
    }

    def 'should report no double' () {
        when:
        GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', zone:'y']])
        then:
        def err =thrown(AbortOperationException)
        and:
        err.message.startsWith('You can\'t specify both \'google.zone\' and \'google.region\' configuration parameters')
    }

}
