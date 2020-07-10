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

import spock.lang.Unroll

import nextflow.exception.AbortOperationException
import nextflow.util.MemoryUnit
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
        !config.preemptible
        !config.disableBinDir
        config.bootDiskSize == null
    }

    def 'should config flags' () {
        when:
        def map = [google:
                           [
                                   project: 'bar',
                                   zone: 'europe-west1-a',
                                   lifeSciences: [preemptible: false, disableRemoteBinDir: false]
                           ]]
        def config = GoogleLifeSciencesConfig.fromSession0(map)
        then:
        config.project == 'bar'
        config.location == 'europe-west2'
        !config.preemptible
        !config.disableBinDir


        when:
        map = [google:
                           [
                                   project: 'bar',
                                   zone: 'europe-west1-a',
                                   lifeSciences: [preemptible: true, disableRemoteBinDir: true]
                           ]]
        config = GoogleLifeSciencesConfig.fromSession0(map)
        then:
        config.project == 'bar'
        config.location == 'europe-west2'
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
        config.location == 'us-central1'
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
        config.location == 'us-central1'
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

    def 'should config boot disk size' () {
        when:
        def config = GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', lifeSciences: [bootDiskSize: "10 GB"]]])
        then:
        config.bootDiskSize == MemoryUnit.of('10 GB')
    }

    def 'should config ssh daemon' () {
        when:
        def config = GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', lifeSciences: [:]]])
        then:
        !config.sshDaemon

        when:
        config = GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', lifeSciences: [sshDaemon:true]]])
        then:
        config.sshDaemon
    }

    def 'should config usePrivateAddress' () {
        when:
        def config = GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', lifeSciences: [:]]])
        then:
        !config.usePrivateAddress

        when:
        config = GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', lifeSciences: [usePrivateAddress:true]]])
        then:
        config.usePrivateAddress
    }

    def 'should config debug mode' () {
        when:
        def config = GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', lifeSciences: [:]]])
        then:
        config.debugMode==null

        when:
        config = GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', lifeSciences: [debug:true]]])
        then:
        config.debugMode == 1

        when:
        config = GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', lifeSciences: [debug:2]]])
        then:
        config.debugMode == 2

    }

    def 'should config copy image' () {
        when:
        def config = GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', lifeSciences: [:]]])
        then:
        config.copyImage == GoogleLifeSciencesConfig.DEFAULT_COPY_IMAGE

        when:
        config = GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', lifeSciences: [copyImage:'foo']]])
        then:
        config.copyImage == 'foo'
    }

    @Unroll
    def 'should return location from region'( ){
        given:
        def config = new GoogleLifeSciencesConfig()

        expect:
        config.bestLocationForRegion(REGION) == LOCATION

        where:
        REGION          | LOCATION
        'europe-west1'  | 'europe-west2'
        'europe-west2'  | 'europe-west2'
        'europe-any'    | 'europe-west2'
        'us-central1'   | 'us-central1'
        'us-any'        | 'us-central1'
        'foo'           | 'us-central1'
    }

    def 'should set requester pays' () {
        when:
        def config = GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', lifeSciences: [:]]])
        then:
        config.enableRequesterPaysBuckets == false

        when:
        config = GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', enableRequesterPaysBuckets:true]])
        then:
        config.enableRequesterPaysBuckets == true

    }
    
    def 'should config cpuPlatform' () {
        when:
        def config = GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', lifeSciences: [:]]])
        then:
        config.cpuPlatform == null

        when:
        config = GoogleLifeSciencesConfig.fromSession0([google:[project:'foo', region:'x', lifeSciences: [cpuPlatform:'Intel Skylake']]])
        then:
        config.cpuPlatform == 'Intel Skylake'

    }

}
