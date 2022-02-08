package nextflow.cloud.azure.config

/**
 * Model Azure azcopy tool config settings from nextflow config file
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */

import groovy.transform.CompileStatic


@CompileStatic
class AzCopyOpts {

    static public final String DEFAULT_BLOCK_SIZE = "4"
    static public final String DEFAULT_BLOB_TIER = "None"

    String blockSize
    String blobTier

    AzCopyOpts() {
        this.blockSize = DEFAULT_BLOCK_SIZE
        this.blobTier =  DEFAULT_BLOB_TIER
    }


    AzCopyOpts(Map config) {
        assert config!=null
        this.blockSize = config.blockSize ?: DEFAULT_BLOCK_SIZE
        this.blobTier = config.blobTier ?: DEFAULT_BLOB_TIER
    }

}
