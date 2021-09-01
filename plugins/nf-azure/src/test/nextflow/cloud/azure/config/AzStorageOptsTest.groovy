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

    def 'should get azure files name, mount point & options'() {
        when:
        def opts1 = new AzStorageOpts([:], [:])

        then:
        opts1.fileName == null
        opts1.relativeMountPoint == null
        opts1.mountOptions == '-o vers=3.0,dir_mode=0777,file_mode=0777,sec=ntlmssp'

        when:
        def opts2 = new AzStorageOpts(
                [fileName: 'source', relativeMountPoint: 'target', 'mountOptions': 'options here'],
                [AZURE_STORAGE_FILE_NAME: 'env-source', AZURE_STORAGE_MOUNT_POINT:'env-target'])
        then:
        opts2.accountName == 'source'
        opts2.accountKey == 'target'
        opts2.mountOptions == 'options here'

        when:
        def opts3 = new AzStorageOpts(
                [:],
                [AZURE_STORAGE_FILE_NAME: 'env-source', AZURE_STORAGE_MOUNT_POINT:'env-target'])
        then:
        opts3.accountName == 'env-source'
        opts3.accountKey == 'env-target'
        opts3.mountOptions == '-o vers=3.0,dir_mode=0777,file_mode=0777,sec=ntlmssp'
    }
}
