/*
 * Copyright 2021, Microsoft Corp
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

package nextflow.cloud.azure.config

import nextflow.Session
import nextflow.util.Duration
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzureConfigTest extends Specification {

    def 'should get azure storage options' () {

        given:
        def KEY = 'xyz1343'
        def NAME = 'container-foo'
        def FILENAME = 'filename-foo'
        def MOUNTPATH = 'mountpath-foo'
        def STORES = 'this-that'
        def SAS = 'foo'
        and:
        def session = Mock(Session) {
            getConfig() >> [ azure:
                                     [storage:[
                                             accountKey: KEY,
                                             accountName: NAME,
                                             fileStores: STORES,
                                             sasToken: SAS,
                                             fileShares: [
                                                     file1: [
                                                             mountPath: '/my/path/1',
                                                             mountOptions:''
                                                     ],
                                                     file2: [
                                                             mountPath: '/my/path/2',
                                                             mountOptions:''
                                                     ]
                                             ]
                                     ] ]]
        }

        when:
        def cfg = AzConfig.getConfig(session)
        then:
        cfg.storage().accountKey == KEY
        cfg.storage().accountName == NAME
        cfg.storage().sasToken == SAS
        cfg.storage().getFileShares().get('file1').mountPath == '/my/path/1'
        cfg.storage().getFileShares().get('file2').mountPath == '/my/path/2'
        and:
        cfg.storage().getEnv() == [AZURE_STORAGE_ACCOUNT_KEY: KEY,
                                   AZURE_STORAGE_ACCOUNT_NAME: NAME,
                                   AZURE_STORAGE_SAS_TOKEN: SAS ]

        and:
        cfg.batch().endpoint == null
        cfg.batch().terminateJobsOnCompletion == true
        cfg.batch().deleteJobsOnCompletion == null
        cfg.batch().deletePoolsOnCompletion == null
        cfg.batch().deleteTasksOnCompletion == null
        cfg.batch().location == null
        cfg.batch().autoPoolMode == null
        cfg.batch().allowPoolCreation == null
        cfg.batch().autoPoolOpts().vmType == 'Standard_D4_v3'
        cfg.batch().autoPoolOpts().vmCount == 1
        cfg.batch().autoPoolOpts().maxVmCount == 3
        cfg.batch().autoPoolOpts().scaleInterval == Duration.of('5 min')
        cfg.batch().autoPoolOpts().autoScale == false
        !cfg.batch().canCreatePool()
    }

    def 'should get azure batch options' () {

        given:
        def KEY = 'xyz1343'
        def NAME = 'container-foo'
        def ENDPOINT = 'http://foo/bar'
        def LOCATION = 'europenorth'
        and:
        def session = Mock(Session) {
            getConfig() >> [ azure:
                                     [batch:[
                                             accountKey: KEY,
                                             accountName: NAME,
                                             endpoint: ENDPOINT,
                                             location: LOCATION,
                                             autoPoolMode: true,
                                             allowPoolCreation: true,
                                             terminateJobsOnCompletion: false,
                                             deleteJobsOnCompletion: true,
                                             deletePoolsOnCompletion: true,
                                             deleteTasksOnCompletion: false,
                                             pools: [ 
                                                myPool: [
                                                    vmType: 'Foo_A1',
                                                    autoScale: true,
                                                    vmCount: 5,
                                                    maxVmCount: 50,
                                                    fileShareRootPath: '/somewhere/over/the/rainbow',
                                                    privileged: true,
                                                    runAs: 'root',
                                                    scaleFormula: 'x + y + z',
                                                    scaleInterval:  '15 min',
                                                    schedulePolicy: 'pack',
                                                    startTask: [ script: 'echo hello-world', privileged: true ]
                                                ]
                                            ]
                                        ]
                                    ]
                                ]
        }

        when:
        def cfg = AzConfig.getConfig(session)
        then:
        cfg.batch().accountKey == KEY
        cfg.batch().accountName == NAME
        cfg.batch().endpoint == ENDPOINT
        cfg.batch().location == LOCATION
        cfg.batch().autoPoolMode == true
        cfg.batch().allowPoolCreation == true
        cfg.batch().terminateJobsOnCompletion == false
        cfg.batch().deleteJobsOnCompletion == true
        cfg.batch().deletePoolsOnCompletion == true
        cfg.batch().deleteTasksOnCompletion == false
        cfg.batch().canCreatePool()
        and:
        cfg.batch().pool('myPool').vmType == 'Foo_A1'
        cfg.batch().pool('myPool').autoScale == true
        cfg.batch().pool('myPool').schedulePolicy == 'pack'
        cfg.batch().pool('myPool').vmCount == 5
        cfg.batch().pool('myPool').maxVmCount == 50
        cfg.batch().pool('myPool').fileShareRootPath == '/somewhere/over/the/rainbow'
        cfg.batch().pool('myPool').scaleFormula == 'x + y + z'
        cfg.batch().pool('myPool').scaleInterval == Duration.of('15 min')
        cfg.batch().pool('myPool').privileged == true
        cfg.batch().pool('myPool').runAs == 'root'
        cfg.batch().pool('myPool').startTask.script == 'echo hello-world'
        cfg.batch().pool('myPool').startTask.privileged == true
    }

    def 'should get azure batch endpoint from account and location' () {

        given:
        def KEY = 'xyz1343'
        def NAME = 'nfbucket'
        def LOCATION = 'europenorth'
        and:
        def session = Mock(Session) {
            getConfig() >> [ azure:
                                     [batch:[
                                             accountKey: KEY,
                                             accountName: NAME,
                                             location: LOCATION ]] ]
        }

        when:
        def cfg = AzConfig.getConfig(session)
        then:
        cfg.batch().accountKey == KEY
        cfg.batch().accountName == NAME
        cfg.batch().endpoint == 'https://nfbucket.europenorth.batch.azure.com'
        cfg.batch().location == LOCATION
    }

    def 'should get azure batch file share root point' () {

        given:
        def KEY = 'xyz1343'
        def NAME = 'container-foo'
        def ENDPOINT = 'http://foo/bar'
        def LOCATION = 'europenorth'

        when:
        def session = Mock(Session) {
            getConfig() >> [ azure:
                                 [batch:[
                                     accountKey: KEY,
                                     accountName: NAME,
                                     endpoint: ENDPOINT,
                                     location: LOCATION,
                                     pools: [
                                         myPool1: [ fileShareRootPath: '/somewhere/over/the/rainbow' ],
                                         myPool2: [ sku:'batch.node.centos 8' ],
                                         myPool3: [ sku:'batch.node.ubuntu 20.04' ],
                                         myPool4: [ sku:'batch.node.debian 10', fileShareRootPath: '/mounting/here' ],
                                         myPool5: [ : ]]
                                 ]] ]
        }
        def cfg = AzConfig.getConfig(session)
        then:
        cfg.batch().pool('myPool1').fileShareRootPath == '/somewhere/over/the/rainbow'
        cfg.batch().pool('myPool2').fileShareRootPath == '/mnt/resource/batch/tasks/fsmounts'
        cfg.batch().pool('myPool3').fileShareRootPath == '/mnt/batch/tasks/fsmounts'
        cfg.batch().pool('myPool4').fileShareRootPath == '/mounting/here'
        cfg.batch().pool('myPool5').fileShareRootPath == AzPoolOpts.DEFAULT_SHARE_ROOT_PATH
   }

    def 'should get azcopy options' () {

        when:
        def session = Mock(Session) {
            getConfig() >> [ azure:
                                     [azcopy:[
                                             blobTier: "Hot",
                                             blockSize: "100" ]] ]
        }

        and:
        def cfg = AzConfig.getConfig(session)

        then:
        cfg.azcopy().blobTier == "Hot"
        cfg.azcopy().blockSize == "100"
    }
}
