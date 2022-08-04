package nextflow.cloud.azure.config

import spock.lang.Specification

/**
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
class AzCopyOptsTest extends Specification {

    def 'should get the default opts'() {
        when:
        def opts1 = new AzCopyOpts([:])
        then:
        opts1.blobTier == 'None'
        opts1.blockSize == '4'
        opts1.overwrite == "false"
        opts1.checkMD5 == "FailIfDifferent"
        opts1.putMD5 == false
    }

    def 'should get the custom opts'() {
        when:
        def opts2 = new AzCopyOpts([blobTier : 'Cool',
                                    blockSize: '100',
                                    putMD5   : true,
                                    overwrite: "true",
                                    checkMD5 : "NoCheck"])
        then:
        opts2.blobTier == 'Cool'
        opts2.blockSize == '100'
        opts2.overwrite == "true"
        opts2.checkMD5 == "NoCheck"
        opts2.putMD5 == true
    }
}
