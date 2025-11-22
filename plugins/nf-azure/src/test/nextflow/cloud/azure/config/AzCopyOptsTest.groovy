package nextflow.cloud.azure.config

import spock.lang.Specification

/**
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
class AzCopyOptsTest extends Specification {

    def 'should get block size and blob tier'() {
        when:
        def opts1 = new AzCopyOpts([:])
        then:
        opts1.blobTier == 'None'
        opts1.blockSize == '4'

        when:
        def opts2 = new AzCopyOpts(
                [blobTier: 'Hot', blockSize: '100'])
        then:
        opts2.blobTier == 'Hot'
        opts2.blockSize == '100'

    }

}
