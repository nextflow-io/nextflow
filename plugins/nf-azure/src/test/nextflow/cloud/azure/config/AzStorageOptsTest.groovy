/*
 * Copyright 2013-2026, Seqera Labs
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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzStorageOptsTest extends Specification {

    def 'should get account name & key'() {
        when:
        def opts1 = new AzStorageOpts([:], [:])
        then:
        opts1.accountName == null
        opts1.accountKey == null

        when:
        def opts2 = new AzStorageOpts(
                [accountName: 'xyz', accountKey: '123'],
                [AZURE_STORAGE_ACCOUNT_NAME: 'env-name', AZURE_STORAGE_ACCOUNT_KEY:'env-key'])
        then:
        opts2.accountName == 'xyz'
        opts2.accountKey == '123'


        when:
        def opts3 = new AzStorageOpts(
                [:],
                [AZURE_STORAGE_ACCOUNT_NAME: 'env-name', AZURE_STORAGE_ACCOUNT_KEY:'env-key'])
        then:
        opts3.accountName == 'env-name'
        opts3.accountKey == 'env-key'
    }

    def 'should get azure files name, root path & options'() {

        given:
        def filesOpts = new LinkedHashMap<String,Map>()
        def storageOpts

        when:
        filesOpts.clear()
        filesOpts['file1'] = [mountOptions: 'mountOptions1', mountPath: 'mountPath1']
        filesOpts['file2'] = [mountOptions: 'mountOptions2', mountPath: 'mountPath2']
        storageOpts = new AzStorageOpts(['fileShares':filesOpts], [:])

        then:
        storageOpts.fileShares.size() == 2
        storageOpts.fileShares.get('file1').getMountPath() == 'mountPath1'
        storageOpts.fileShares.get('file1').getMountOptions() == 'mountOptions1'
        storageOpts.fileShares.get('file2').getMountPath() == 'mountPath2'
        storageOpts.fileShares.get('file2').getMountOptions() == 'mountOptions2'

        when:
        filesOpts.clear()
        filesOpts['file1'] = [mountPath: 'mountPath1']
        storageOpts = new AzStorageOpts(['fileShares':filesOpts], [:])

        then:
        storageOpts.fileShares.size() == 1
        storageOpts.fileShares.get('file1').getMountPath() == 'mountPath1'
        storageOpts.fileShares.get('file1').getMountOptions() == AzFileShareOpts.DEFAULT_MOUNT_OPTIONS
    }
}
